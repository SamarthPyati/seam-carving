#include <assert.h>

#define NOB_STRIP_PREFIX
#define NOB_IMPLEMENTATION
#include "external/nob.h"

#include <time.h>

// #define NO_OPTIMISATION
#define COMPILE_NATIVE

static double get_time(void)
{
    struct timespec tp = {0};
    int ret = clock_gettime(CLOCK_MONOTONIC, &tp);
    assert(ret == 0);
    return tp.tv_sec + tp.tv_nsec * 1e-9;
}

static void incorporate_omp(Cmd *cmd) {
#if defined(__clang__)
    cmd_append(cmd, "-Xpreprocessor", "-fopenmp");
#else 
    cmd_append(cmd, "-fopenmp");
#endif
    cmd_append(cmd, "-I/opt/homebrew/opt/libomp/include");
    cmd_append(cmd, "-L/opt/homebrew/opt/libomp/lib");
    cmd_append(cmd, "-lomp");
}

static void cc(Cmd *cmd)
{
    cmd_append(cmd, "cc");
    cmd_append(cmd, "-Wall", "-Wextra", "-Wno-unused-function", "-g", "-std=c17");
    cmd_append(cmd, "-fsanitize=address,undefined");

    incorporate_omp(cmd);

#if defined(COMPILE_NATIVE)
    cmd_append(cmd, "-march=native");
#endif

#if !defined(NO_OPTIMISATION)
    cmd_append(cmd, "-O3");
#else
    cmd_append(cmd, "-O0");
#endif
}

static bool build_binary(Cmd *cmd, const char *source, const char *output)
{
    if (needs_rebuild1(output, source)) {
        cmd->count = 0;
        cc(cmd);
        cmd_append(cmd, "-I./include");
        cmd_append(cmd, "-I./external");
        cmd_append(cmd, "-o", output);
        cmd_append(cmd, source);
        cmd_append(cmd, "-lm");
        return cmd_run(cmd);
    }

    nob_log(INFO, "%s is up to date", output);
    return true;
}

int main(int argc, char **argv)
{
    GO_REBUILD_URSELF(argc, argv);

    const char *program = shift_args(&argc, &argv);
    (void)program;

    Cmd cmd = {0};

    if (!mkdir_if_not_exists("./build/")) return 1;

    const char *main_source = "src/seam_carving.c";
    const char *main_output = "./build/seam-carving";
    if (!build_binary(&cmd, main_source, main_output)) return 1;

    cmd.count = 0;
    cmd_append(&cmd, main_output);
    da_append_many(&cmd, argv, argc);
    double begin = get_time();
    if (!cmd_run(&cmd)) return 1;
    double end = get_time();
    nob_log(INFO, "Execution time: %lf second(s)", end - begin);

    return 0;
}