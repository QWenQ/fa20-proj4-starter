from distutils.core import setup, Extension
import sysconfig

def main():
    CFLAGS = ['-g', '-Wall', '-std=c99', '-fopenmp', '-mavx', '-mfma', '-pthread', '-O3']
    LDFLAGS = ['-fopenmp']
    # Use the setup function we imported and set up the modules.
    # You may find this reference helpful: https://docs.python.org/3.6/extending/building.html
    # TODO: YOUR CODE HERE
    module1 = Extension(name = 'numc',
                        sources = ['numc.c'])
    
    setup(name = 'numc',
          description = 'DIY numc module',
          ext_modules = [module1])

if __name__ == "__main__":
    main()
