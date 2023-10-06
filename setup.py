from setuptools import find_packages, setup
from setuptools.command.install import install as _install
from setuptools.command.develop import develop as _develop

import importlib
import logging
import subprocess
import distutils.cmd


class DepsCommand(distutils.cmd.Command):
    """Custom command to install dependencies only, this
    is intended to be used to take advantage of Docker's
    caching.
    This installs into the user home so don't run it locally, use
    "pip install .[dev]" inside a virtual environment instead
    """

    description = "Install dependencies"
    user_options = [
        # The format is (long option, short option, description).
        ("env=", None, "Name of environment")
    ]

    def initialize_options(self):
        self.env = ""

    def finalize_options(self):
        if self.env:
            assert self.env in ["dev"], "env: {} not available".format(self.env)

    def run(self):
        command = ["pip", "install", "--user"]
        command.extend(install_requires)
        if self.env:
            command.extend(extras["dev"])
        self.announce(
            "Running command: %s" % " ".join(command), level=distutils.log.INFO
        )
        subprocess.check_call(command)


def _safe_read_lines(f):
    with open(f) as in_f:
        r = in_f.readlines()
    r = [l.strip() for l in r]
    return r


console_scripts = [
    "simulate-cancer-cells=neoag_dt.cli.simulate_neoag_cancer_cells:main",
    "optimize-vaccine-ilp=neoag_dt.cli.optimize_vaccine_ilp:main",
    "create-bar-chart=neoag_dt.cli.create_selected_vaccine_element_bar_chart:main",
    "evaluate-vaccine-response=neoag_dt.cli.evaluate_vaccine_response:main"
]

install_requires = _safe_read_lines("./requirements.txt")

tests_require = [
    'pytest',
    'coverage',
    'pytest-cov',
    'coveralls',
    'pytest-runner',
    'pylint'
]

gpu_requires = []

docs_require = [
    'sphinx==3.5.4',
    'sphinx_rtd_theme',
    'Jinja2<3.1'
]

all_requires = (
    tests_require +
    gpu_requires +
    docs_require
)

extras = {
    'test': tests_require,
    'gpu': gpu_requires,
    'docs': docs_require,
    'all': all_requires
}

classifiers = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Programming Language :: Python :: 3 :: Only',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
]


def _post_install(self):
    import site
    importlib.reload(site)


class my_install(_install):
    def run(self):
        level = logging.getLevelName("INFO")
        logging.basicConfig(level=level,
            format='%(levelname)-8s : %(message)s')

        _install.run(self)
        _post_install(self)


class my_develop(_develop):  
    def run(self):
        level = logging.getLevelName("INFO")
        logging.basicConfig(level=level,
            format='%(levelname)-8s : %(message)s')

        _develop.run(self)
        _post_install(self)


def readme():
    with open('README.md') as f:
        return f.read()


def description():
    description = ("This package contains a digital twin approach for "
        "designing personalized cancer vaccines.")
    return description


setup(
    name='neoag-dt',
    version='1.0',
    description=description(),
    long_description=readme(),
    keywords="peptide epitope immunogenicity mhc hla ",
    author="Filippo Grazioli, Anja Moesch - NEC Labs Europe",
    author_email="anja.moesch@neclab.eu",
    license='BSD 3-clause "New" or "Revised License"',
    packages=find_packages(),
    install_requires=install_requires,
    cmdclass={
        'install': my_install,  # override install
        'develop': my_develop,   # develop is used for pip install -e .
        'deps': DepsCommand  # for Docker
    },
    include_package_data=True,
    tests_require=tests_require,
    extras_require=extras,
    entry_points={
        'console_scripts': console_scripts
    },
    zip_safe=False,
    classifiers=classifiers,    
)
