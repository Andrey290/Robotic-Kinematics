cmake_minimum_required(VERSION 3.16)
project(RobotControl
        VERSION 0.1
        LANGUAGES C)

# -----------------------
# Настройки проекта
# -----------------------
# Требуем стандарт C11
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Рекомендуется включать экспорт compile_commands для редакторов (clangd)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Режим сборки по умолчанию (для single-config генераторов)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type" FORCE)
endif()

# -----------------------
# Файлы проекта
# -----------------------
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(INC_DIR ${CMAKE_SOURCE_DIR}/inc)

set(SOURCES
    ${SRC_DIR}/main.c
    ${SRC_DIR}/kinematic.c
)

set(HEADERS
    ${INC_DIR}/globalconstants.h
    ${INC_DIR}/kinematic.h
    ${INC_DIR}/mathutils.h
)

# -----------------------
# Цели
# -----------------------
add_executable(robotcontrol ${SOURCES} ${HEADERS})

if(UNIX)
  target_link_libraries(robotcontrol PRIVATE m)
endif()

# Include dir — используем target_* чтобы зависимость была локальной для цели
target_include_directories(robotcontrol
    PRIVATE
        ${INC_DIR}
)

# Warnings — кроссплатформенно
if(MSVC)
    target_compile_options(robotcontrol PRIVATE /W4 /permissive-)
else()
    target_compile_options(robotcontrol PRIVATE -Wall -Wextra -Wpedantic)
    # для разработки можно добавить -Werror (осторожно)
    # target_compile_options(robotcontrol PRIVATE -Werror)
endif()

# Опционально: defines
target_compile_definitions(robotcontrol PRIVATE DEFAULT_N_LINKS=6 DH_COLS=4 _USE_MATH_DEFINES)

# Удобные цели для установки/тестов можно добавить позже

