#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define STMT(S) do { S } while (0)

#define LOG_ERROR(fmt, ...) STMT(                                                   \
    fprintf(stderr, "%s:%d: ERROR: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__);  \
)

#define LOG_ERROR_AND_ABORT(fmt, ...) STMT(                                         \
    fprintf(stderr, "%s:%d: ERROR: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__);  \
    abort();                                                                        \
)

#define LOG_FATAL(fmt, ...) STMT(                                                   \
    fprintf(stderr, "%s:%d: FATAL: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__);  \
    abort();                                                                        \
)

/* Usage: LOG_USAGE(argv[0], "<image_file_path> <seams_to_remove>"); or with format args */
#define LOG_USAGE(prog, fmt, ...) STMT(                                 \
    fprintf(stderr, "USAGE: %s " fmt "\n", (prog), ##__VA_ARGS__);      \
    exit(EXIT_FAILURE);                                                 \
)


#define CLAMP(val, lo, hi) ( ((val) > (hi)) ? (hi) : (((val) < (lo)) ? (lo) : (val)) )

/* in-place clamp that avoids double-evaluation (GCC/Clang typeof), fallback is simple assign */
#ifdef __GNUC__
#define CLAMP_ASSIGN(v, lo, hi) STMT(                       \
    __typeof__(v) _v = (v);                                 \
    __typeof__(lo) _lo = (lo);                              \
    __typeof__(hi) _hi = (hi);                              \
    (v) = (_v > _hi) ? _hi : ((_v < _lo) ? _lo : _v);       \
)
#else
#define CLAMP_ASSIGN(v, lo, hi) STMT( (v) = CLAMP(v, lo, hi); )
#endif

#endif // UTILS_H_