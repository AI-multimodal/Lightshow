from setuptools import setup

setup(
    name='xanes_bench',
    version='0.01',
    packages=['runScripts'],
    url='',
    license='MIT License',
    author='xiaqu',
    author_email='john.vinson@nist.gov; dlu@bnl.gov; vorwerk@physik.hu-berlin.de; xiaqu@bnl.gov',
    description='Benchmark code and data for XANES Simulation',
    entry_points={
        "console_scripts": [
            "fetchSingle = runScripts.fetchSingle:main",
            "makeXspectraInputs = runScripts.makeXspectraInputs"
        ]
    }
)
