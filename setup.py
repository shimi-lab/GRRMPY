from setuptools import setup, find_packages

console_scripts = [
    'read_irc=grrmpy.command.read_irc:read_irc',
    'read_lup=grrmpy.command.read_lup:read_lup',
]

setup(
    name='grrmpy',
    version='1.0',
    packages=find_packages(),
    entry_points={'console_scripts': console_scripts},
)
