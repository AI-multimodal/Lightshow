from setuptools import setup


with open("requirements.txt") as f_req:
    required_list = [line.rstrip() for line in f_req.readlines()]


setup(
    name='xanes_bench',
    version='0.01',
    packages=['runScripts'],
    url='',
    license='MIT License',
    author='xiaqu',
    author_email='john.vinson@nist.gov; dlu@bnl.gov; vorwerk@physik.hu-berlin.de; xiaqu@bnl.gov',
    description='Benchmark code and data for XANES Simulation',
    python_requires='>=3.7',
    install_requires=required_list,
    entry_points={
        "console_scripts": [
            "fetchSingle = runScripts.fetchSingle:main",
            "makeXspectraInputs = runScripts.makeXspectraInputs"
        ]
    }
)