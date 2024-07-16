//
// Created by cepheid on 5/7/24.
//

#ifndef GAME_ENGINE_ASSERTS_H
#define GAME_ENGINE_ASSERTS_H

#ifdef DEBUG

#include <cstdio>
// #include "dbg.h"

#define ANSI_RED     "\x1b[31m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_YELLOW  "\x1b[33m"
#define ANSI_BLUE    "\x1b[34m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_CYAN    "\x1b[36m"
#define ANSI_RESET   "\x1b[0m"

// inline std::FILE* LOG_FILE_HANDLE = std::fopen("logfile.log", "a+");
#define DEBUG_WARN(fmt, ...) std::fprintf(stdout, ANSI_YELLOW "[WARNING] " ANSI_RESET fmt, ##__VA_ARGS__)
#define DEBUG_ERR(fmt, ...) std::fprintf(stderr, fmt, ##__VA_ARGS__)
#define DEBUG_MSG(fmt, ...) std::fprintf(stdout, ANSI_BLUE fmt ANSI_RESET, ##__VA_ARGS__)

#define BREAKPOINT() __asm__ __volatile__("int 3")

/*
// check the expression and fail if it is false
#define ASSERT(expr, fmt, ...) \
if (expr) {} \
else { \
std::fprintf(stderr, "[ASSERT](%s:$d) " #expr " : " fmt, __FILE__, __LINE__, __VA_ARGS__); \
BREAKPOINT(); \
} //static_assert(true) // for that sweet trailing semicolon action
*/

#else

#define ASSERT(expr) // evaluates to nothing
#define SLOW_ASSERT(expr) // Not sure how this is supposed to work

#define LOG_WARNING(fmt, ...)
#define LOG_ERROR(fmt, ...)
#define DEBUG_MSG(fmt, ...)

#endif // DEBUG
#endif //GAME_ENGINE_ASSERTS_H
