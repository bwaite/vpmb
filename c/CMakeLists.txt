project(vpmb)

set (SOURCES
  main.c
  vpmb.c
  cJSON.c
  strlcpy.c
  utils.c
)
set(CMAKE_BUILD_TYPE debug)
#add "-ansi" back before the final release
set(CMAKE_C_FLAGS_DEBUG "-W -Wall -Wextra -Wpedantic -ggdb3 -lm ")
set(CMAKE_C_FLAGS_RELEASE "-O3 -W -Wall -Wextra -Wpedantic -lm ")

add_executable(vpmb ${SOURCES})