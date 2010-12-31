project(vpmb)

set (SOURCES
  main.c
  vpmb.c
  cJSON.c
  strlcpy.c
)
set(CMAKE_BUILD_TYPE debug)
#add "-ansi" back before the final release
set(CMAKE_C_FLAGS_DEBUG "-W -Wall -Wextra -ggdb3 -lm -std=c99")
set(CMAKE_C_FLAGS_RELEASE "-O3 -W -Wall -Wextra -lm -std=c99")

add_executable(vpmb_c ${SOURCES})
