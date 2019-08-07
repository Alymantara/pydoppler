import setuptools

# upload to pip
# pip install .
# python3 setup.py sdist bdist_wheel
# twine upload dist/pydoppler-0.1.8.tar.gz

import os

paths = []
def package_files(directory):
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))

package_files('pydoppler/fortran_code')
package_files('pydoppler/test_data')


setuptools.setup(
     name='pydoppler',
     version='0.2.0',
     packages=['pydoppler'] ,
     package_data={'': paths},
     author="Juan V. Hernandez Santisteban",
     author_email="jvhs1@st-andrews.ac.uk",
     description="A python wrapper for Henk Spruit's doppler tomography sowftare",
   long_description_content_type="text/markdown",
     url="https://github.com/alymantara/pydoppler",
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
         ],
 )
