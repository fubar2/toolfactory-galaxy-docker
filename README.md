# toolfactory-galaxy-docker
Docker Galaxy container with the ToolFactory and demonstration tools

This repository contains the Dockerfile and files needed to build the container image at
quay.io/fubar2/toolfactory-galaxy-docker for building Galaxy tool wrappers inside Galaxy.

**Audience**

Many developers will be happy to use Planemo on the command line, but for users who
prefer building and testing new tools inside the Galaxy UI, this docker container may be useful.
Although Galaxy as an IDE may seem odd, it works pretty well.

**Building**

Using the quay.io image will be simplest, but building an image from this Dockerfile is straightforward.
It is based on https://github.com/bgruening/docker-galaxy-stable so has release 20.09 at present. Run the
supplied workflow using the history data files in the supplied history to generate all the sample tools -this
has already been done in the supplied image.

**Starting the Docker image**

There are two supplied shell scripts with appropriate command lines.
1. runclean.sh will clear out the export/ directory before running.
2. rundirty will not clear out export/

Both will open ports including 8080/9009.
Privileged mode is needed to run the planemo biocontainer as a docker-in-docker image.

**Dependency management**

Conda is used for all dependency management although tools that use system utilities like sed, bash or awk
do not require any dependencies if the utilities are available on job execution nodes. Sed and friends are available as Conda (conda-forge) dependencies if necessary.
The ToolFactory relies on galaxyxml to generate tool xml, and uses ephemeris and
bioblend to load tools to the toolshed and to Galaxy.

**Building new Galaxy tools**

Clearly, it is best to have done all the required planning
to wrap a new script or binary. At the very least, the user will need to know what the command line for the script or
binary should look like, and have some sample inputs in the working history together with values for every required
parameter for these sample inputs, since these must be supplied on the ToolFactory form for preparing
the inbuilt tool test.

A new tool is specified by filling in the usual Galaxy tool form. Small sample data and test parameter settings must be
supplied. Once the form is completed, executing the ToolFactory
will build a new XML tool wrapper including a functional test based on the sample settings and data.

If the Planemo test passes, the tool can be optionally uploaded to the local Galaxy used in the image for more testing.

A local toolshed runs inside the container to allow an automated installation, although any toolshed can be specified
for this process.

Steps in building a new Galaxy tool are all conducted through Galaxy running in the docker container:

1. Open the ToolFactory tool in the Galaxy running in the container at http://localhost:8080
2. Fill in the form - supply the new tool name, tool model, dependencies, i/o and sample files, parameters and sample settings, optional script.
3. Execute the tool to create a new XML tool wrapper using the sample inputs and parameter settings for the inbuilt tool test. Planemo is run to generate the outputs
    from the test. The complete toolshed archive is written to the history together with the planemo test report. Optionally the new tool archive can be uploaded
    to the toolshed running in the same container (http://localhost:9009) and then installed inside the Galaxy in the container for further testing.
4. If the test fails, rerunning the failed history job allows errors on the tool form to be edited before rerunning until everything works correctly.

**Outputs**

New tools are ordinary Galaxy tools, ready to publish in any Galaxy Toolshed and install to any running Galaxy instance.
They are fully workflow compatible and work exactly like any hand-written tool.

**Models for tool execution**

The simplest tool model wraps a simple script or Conda dependency package requiring zero parameters, such as filters that take input from STDIN and write to STDOUT.
These can be configured to take STDOUT and write it to a new history item, and to read a user selected history data file on STDIN.

This simple model can be explored in the Tacrev demonstration tool found in the Galaxy running in the container. It takes an input file, uses a bash script to
read that file on STDIN, run the unix tac utility to reverse catenate the lines, and pipe that output to the unix rev utility to reverse each text line.

That simple model can be made more complicated by adding more input data files, parameters and more output files, to allow more complicated
scripts or dependent binaries that require:

1. An arbitrary number of input data files selected by the user from existing history data
2. An arbitrary number of output data files written to the user's history
3. An arbitrary number of user supplied parameters. These can be passed as command line arguments to the script or the dependency package. Either
positional or named (argparse) style command line parameter passing can be used.

More complex models can be seen in the Sedtest, Pyrevpos and Pyrevargparse tools illustrating positional and argparse parameter passing.

The most complex demonstration is the Planemo advanced tool tutorial BWA tool. There is one version using a command-override to implement
exactly the same command structure in the Planemo tutorial. A second version uses a bash script and positional parameters to achieve the same
result. Some users may find the bash version more familiar and cleaner but the choice is yours.

**Requirements**

These are all managed automagically. Planemo is used for testing and runs in a biocontainer currently at quay.io/fubar2/planemo-biocontainer.
This is needed because at present, Planemo seems to have a bug allowing it to leak dependencies back into the calling environment leaving that
environment permanently damaged.  So, it is run completely isolated in a separate container. The docker python SDK is used to manage the
complexities of running docker-in-docker inside the running ToolFactory tool. Trust me - there are complications.

**Caveats**

This docker image requires privileged mode so exposes all sorts of security risks. Please, do not run it in any situation where that is
a problem - never, ever on a public facing Galaxy server. On a laptop or workstation should be fine in a non-hostile environment.









