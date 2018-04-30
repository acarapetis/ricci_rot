def FlagsForFile( filename, **kwargs ):
    #if re.match(r'.*\.[ch]', filename) is not None:
    return {
        'flags': [ '-x', 'c', '-Werror', '-lglut', '-lGL', '-lGLU', '-lm', '-isystem', '/usr/include', '-isystem', '/usr/include/GL',
            '-isystem', '/home/anthony/src/emsdk/emscripten/1.37.38/system/include'
        ],
    }
