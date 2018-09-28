from setuptools import setup

setup(name='graphalign',
      version='0.1',
      description='Graph Align',
      url='http://github.com/uio-bmi/graphalign',
      author='Ivar Grytten and Knut Rand',
      author_email='knutdrand@gmail.com',
      license='MIT',
      packages=['graphalign'],
      zip_safe=False,
      install_requires=['numpy', 'pytest'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points={
        'console_scripts': ['graphalign=graphalign.interface:main']}
      )
