#define NOB_IMPLEMENTATION
#include "external/nob.h"


#define COMPILE_FASTER 

int main(int argc, char *argv[])
{   
    NOB_GO_REBUILD_URSELF(argc, argv);
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc");
    nob_cmd_append(&cmd, "-Wall", "-Wextra", "-ggdb", "-std=c17");

#if defined(COMPILE_FASTER)
    nob_cmd_append(&cmd, "-O3");
#else 
    nob_cmd_append(&cmd, "-O0");
#endif
    nob_cmd_append(&cmd, "-o", "main");
    nob_cmd_append(&cmd, "stb_image.o", "stb_image_write.o", "nob.o", "main.c");

    if (!nob_cmd_run_sync(cmd)) return -1;

    // cmd.count = 0;
    // nob_cmd_append(&cmd, "./main");
    // if (!nob_cmd_run_sync(cmd)) return -1;
    return 0;
}