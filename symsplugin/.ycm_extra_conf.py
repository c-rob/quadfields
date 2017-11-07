# Minimal compilation flags gcc-like for dsPIC C project

def FlagsForFile( filename, **kwargs ):
	
	localIncludePaths = [
		'.',
		'/home/roberto/bin/V-REP/programming/include/'
	]
	defines = [
		'__linux'
	]

	baseCFlags = ['-x', 'c++', '-std=c++11', '-Wall', '-Wextra', '-Werror',
			'-g', '-O0']
	includeFlags = ['-I'+p for p in localIncludePaths]
	definesFlags = ['-D'+p for p in defines]

	return {
	'flags': baseCFlags + includeFlags + definesFlags
	}
