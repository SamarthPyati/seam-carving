#include <assert.h>

#define NOB_IMPLEMENTATION
#include "external/nob.h"

#include <time.h>

// #define COMPILE_FASTER
#define COMPILE_NATIVE

static double get_time(void)
{
    struct timespec tp = {0};
    int ret = clock_gettime(CLOCK_MONOTONIC, &tp);
    assert(ret == 0);
    return tp.tv_sec + tp.tv_nsec * 1e-9;
}

static void cc(Nob_Cmd *cmd)
{
    nob_cmd_append(cmd, "cc");
    nob_cmd_append(cmd, "-Wall", "-Wextra", "-Wno-unused-function", "-g", "-std=c17");
    nob_cmd_append(cmd, "-fsanitize=address,undefined");

#if defined(COMPILE_NATIVE)
    nob_cmd_append(cmd, "-march=native");
#endif

#if !defined(COMPILE_FASTER)
    nob_cmd_append(cmd, "-O3");
#else
    nob_cmd_append(cmd, "-O0");
#endif
}

static bool build_binary(Nob_Cmd *cmd, const char *source, const char *output)
{
    if (nob_needs_rebuild1(output, source)) {
        cmd->count = 0;
        cc(cmd);
        nob_cmd_append(cmd, "-I./include");
        nob_cmd_append(cmd, "-I./external");
        nob_cmd_append(cmd, "-o", output);
        nob_cmd_append(cmd, source);
        nob_cmd_append(cmd, "-lm");
        return nob_cmd_run_sync(*cmd);
    }

    nob_log(NOB_INFO, "%s is up to date", output);
    return true;
}

int main(int argc, char **argv)
{
    NOB_GO_REBUILD_URSELF(argc, argv);

    const char *program = nob_shift_args(&argc, &argv);
    (void)program;

    Nob_Cmd cmd = {0};

    if (!nob_mkdir_if_not_exists("./build/")) return 1;

    const char *main_source = "src/main.c";
    const char *main_output = "./build/seam-carving";
    if (!build_binary(&cmd, main_source, main_output)) return 1;

    cmd.count = 0;
    nob_cmd_append(&cmd, main_output);
    nob_da_append_many(&cmd, argv, argc);
    double begin = get_time();
    if (!nob_cmd_run_sync(cmd)) return 1;
    double end = get_time();
    nob_log(NOB_INFO, "Execution time: %lf second(s)", end - begin);

    return 0;
}