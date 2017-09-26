# hts-db
The hts-db is a simple but powerful system for storing and analyzing chemical-bioactivity assay data. The backend is a mongodb NoSQL database with a python mongoengine object-oriented interface. All data analysis and visualization is handled in python using Pandas, numpy, scipy and matplotlib. 

# Installation

## Requirements

### Operating system 
This system has only been tested under [Red Hat Enterprise Linux v6](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) and [Ubuntu Server 16.04 LTS](http://www.ubuntu.com). The Ubuntu system can be downloaded freely and installed on any Intel PC. 

### MongoDB

Install the latest version of [mongodb](http://www.mongodb.com). MongoDB is a document-oriented database that is quite handy for storing complex data structures [here's](https://docs.mongodb.com/getting-started/shell/introduction/) an introduction to MongoDB. The hts-db system uses a Python object-oriented interface to MongoDB called (mongoengine)[http://mongoengine.org/], which really takes the pain out of working with complex persistent data as you'll see in some of the python notebooks in this project. 
 
### Python

This code has been tested on [Python 2.7](http://python.org). The easiest way to install Python 2.7 is via [Anaconda](https://www.continuum.io/downloads). Install the python packages given in py_requirements.txt as follows:

```
unix> pip install -r py_requirements.txt
```

### Jupyter
[Jupyter notebook](http://jupyter.org/) is an interactive computing environment based and is freely available. It can be installed using Anaconda. Please the the documents about Jupyter or find a tutorial on YouTube.  After jupyter is installed read the [quickstart instructions](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/) to create a configuration file. Set the following variables in your jupyter_notebook_config.py:

```python
c.NotebookApp.notebook_dir = your_notebook_directory
c.NotebookApp.port = 7777
```

Start the notebook server from the command line using (I always run this in a [screen](https://www.gnu.org/software/screen/manual/screen.html)):

```
unix> jupyter notebook
```

### Python source: Place the python source files in your PYTHONPATH. 

# Loading MongoDB data
The data required for this analysis is available as a mongodump file and it must be downloaded from [here](https://tinyurl.com/y8skjbds).  Untar this file using the following commmand (will need bunzip2 decompression). Use the [mongorestore](https://docs.mongodb.com/manual/reference/program/mongorestore/) tool (installed as part of the mongodb package) to recreate the htsdb database (by user= user with password = passwd):-

```
unix> mongorestore -u user -p passwd -d htsdb --gzip mongodump 
```


# Testing the system
After you have completed the above requirements save the hts-db-testing.ipynb notebook to your_notebook_directory and open [this page](http://localhost:7777) in your browser. Each section of code (called a cell) can be executed by pressing the ">|" button on the menu or by the ctrl-enter key combination to produce the output. All cells must be run sequentially in order to ensure the variables have been initialized. If you can generate the concentration response heatmaps then the above process was successful .... 


