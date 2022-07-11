from setuptools import setup

setup(
    name='dCas9_fusions',
    version='0.0.1',

    author='Jeff Hussmann',
    author_email='jeff.hussmann@gmail.com',

    packages=[
        'dCas9_fusions',
    ],

    python_requires='>=3.7',

    install_requires=[
        'tqdm',
        'knock_knock>=0.3.8',
    ],

)
