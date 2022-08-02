from setuptools import setup
from pathlib import Path

targets_fns = []
targets_dir = Path('dCas9_fusions/targets')
for fn in targets_dir.glob('**/*'):
    targets_fns.append(str(fn.relative_to('dCas9_fusions')))

setup(
    name='dCas9_fusions',
    version='0.0.1',

    author='Jeff Hussmann',
    author_email='jeff.hussmann@gmail.com',

    packages=[
        'dCas9_fusions',
    ],

    package_data={
        'dCas9_fusions': targets_fns,
    },

    scripts=[
        'dCas9_fusions/count_dCas9_fusion_domains',
    ],

    python_requires='>=3.7',

    install_requires=[
        'tqdm',
        'hits>=0.3.3',
        'knock_knock>=0.3.8',
    ],

)
