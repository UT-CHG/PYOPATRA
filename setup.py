# Largely taken from https://www.benjack.io/2018/02/02/python-cpp-revisited.html and https://github.com/pybind/cmake_example

import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools.command.develop import develop


# From https://stackoverflow.com/questions/18725137/how-to-obtain-arguments-passed-to-setup-py-from-pip-with-install-option
class CommandMixin:
    user_options = install.user_options + [
        ('cxxcompiler=', None, None),
        ('sandybridge', None, 'Compile for Sandy Bridge Architectures')
    ]

    def initialize_options(self):
        super().initialize_options()
        self.cxxcompiler = None
        self.sandybridge = None

    def finalize_options(self):
        print("value of cxxcompiler is", self.cxxcompiler)
        print("use sandy bridge is ", self.sandybridge)
        super().finalize_options()

    def run(self):
        global cxxcompiler
        global sandybridge
        cxxcompiler = self.cxxcompiler # will be 1 or None
        sandybridge = self.sandybridge

        super().run()

class InstallCommand(CommandMixin, install):
    user_options = getattr(install, 'user_options', []) + CommandMixin.user_options

class DevelopCommand(CommandMixin, develop):
    user_options = getattr(develop, 'user_options', []) + CommandMixin.user_options

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                                   out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DBUILD_PYTHON=ON']

        try:
            if cxxcompiler:
                cmake_args.append('-DCMAKE_CXX_COMPILER=' + cxxcompiler)
        except NameError:
            pass

        try:
            if sandybridge:
                cmake_args.append('-DSANDY_BRIDGE=ON')
        except NameError:
            pass

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)
        print()  # Add an empty line for cleaner output

setup(
    name='PYOPATRA',
    version='0.1',
    author='Georgia Stuart',
    author_email='georgia.stuart@austin.utexas.edu',
    description='A lagrangian particle tracker (LPT) open ocean solver',
    long_description='',
    cmdclass={
        'install': InstallCommand,
        'develop': DevelopCommand,
        'build_ext': CMakeBuild
    },
    # tell setuptools to look for any packages under 'src'
    packages=find_packages('src'),
    # tell setuptools that all packages will be under the 'src' directory
    # and nowhere else
    package_dir={'':'src'},
    # add an extension module named 'python_cpp_example' to the package
    # 'python_cpp_example'
    ext_modules=[CMakeExtension('PYOPATRA/PYOPATRA')],
    # add custom build_ext command
    zip_safe=False,
)