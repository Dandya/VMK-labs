cmake_minimum_required(VERSION 3.12)  # Минимальная версия CMake
project(labs_vmk)                     # Название проекта

set(CMAKE_CXX_STANDARD 17)            # Используем C++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)   # Требуем поддержку C++17

# Включаем OpenMP (если доступно)
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Оптимизация (-O2)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O2)
endif()

# Исходные файлы
set(SOURCES
    main.cpp
    labs.cpp
    lab1.cpp
    lab2.cpp
    lab3.cpp
    # Lab4.cpp
)

# Создаем исполняемый файл
add_executable(labs ${SOURCES})

# Подключаем OpenMP (если найден)
if(OpenMP_FOUND)
    target_link_libraries(labs PRIVATE OpenMP::OpenMP_CXX)
endif()