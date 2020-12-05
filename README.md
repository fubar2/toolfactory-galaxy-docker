# toolfactory-galaxy-docker

## Docker Galaxy container with the ToolFactory and demonstration tools.
### A Galaxy IDE for tool wrappers.

## Overview

This repository contains the Dockerfile and files needed to build and run the container image available from
quay.io/fubar2/toolfactory-galaxy-docker. It allows the user to build new Galaxy tool wrappers inside Galaxy,
using the specialised ToolFactory tool - a Galaxy tool that generates tool wrappers.

It relies on https://github.com/bgruening/docker-galaxy-stable (20.09 at present) and
the ToolFactory from https://toolshed.g2.bx.psu.edu/view/fubar/tool_factory_2. Testing is done
using Planemo from https://github.com/galaxyproject/planemo

## Security warning for this container.

*This container must be run in privileged mode, exposing important potential security risks*

The ToolFactory is a specialised tool that requires a privileged docker container, based on docker-galaxy-stable. Although
the ToolFactory tool will only run for administrative users, you are strongly advised not to expose this container on any
public facing production hardware because that potential opportunity for privilege escalation is generally not an acceptable risk.

All generated tools are just normal Galaxy tool XML wrappers and tests, with no additional security vulnerabilities from the ToolFactory itself.

The ToolFactory just makes writing tools easier and uses a familiar Galaxy tool interface. Like any Galaxy tool, ToolFactory products *could* be
constructed with malicious code. Outside the privileged docker container used to generate them, generated tools run as normal Galaxy jobs
in an isolated environment where damage will be limited.

## Intended audience and users

Most developers use Planemo on the command line to prepare and test new Galaxy tools. This image
is for users who prefer building and testing new tools *inside the Galaxy UI* !

## Building the container

Building an image from this Dockerfile is straightforward but most users will probably want to run one of the startup
scripts supplied. They will download the https://quay.io/fubar2/toolfactory_galaxy_docker image or build it, depending
on your preference.

Running the image from a directory with a subdirectory (export/) is recommended for Bjoern's Galaxy docker because
it allows persistence between container restarts.

Occasionally, Galaxy fails to launch because the database is not available and after a minute, I just stop that and
start afresh. Usually launches fine then. The script does a thorough cleanup before it starts. Maybe the docker
daemon needs more frequent restarts.

There are two supplied shell scripts with appropriate command lines to either start a fresh clean instance or restart
the image without clearing out the old database and tools:

1. runclean.sh will clear out the export/ directory before running. The first run may take 10-15 minutes to configure and start
Galaxy, a toolshed and install the ToolFactory tool with dependencies.

2. rundirty will *not* empty the export/ directory, so can start a previous instance although using docker stop and start is the best way
to preserve a working container.

Both will open ports including 8080/9009. Privileged mode is needed to run the planemo biocontainer as a docker-in-docker image.
If this is the first time you build the image, please immediately run the restart.sh script. The workflow will break otherwise for reasons
unknown. A restart seems to work.

`docker exec [container id] sh /usr/local/bin/restartall.sh`

Wait until activity dies down before logging in.

## Demonstration tools

There are currently 8 tools plus the ToolFactory that are built by running a workflow included in the image.
They illustrate the range of models for tool execution that the ToolFactory can produce, described in the next section.
Run the workflow, selecting supplied history dataset names as shown to match the filename prompts in each input file field.
Why doesn't the workflow runner do that?:

![Input file configuration to run the supplied workflow](files/TFWorkflow_setup.png?raw=true "Input file configuration to run the supplied workflow")

The history in which they were all built allows each one to be rerun so you can try re-running them to generating new tools by editing their original build tool form.

See below for generated examples.

## ToolFactory generated tools are ordinary Galaxy tools

A ToolFactory generated tool that passes the Planemo test is ready to publish in any Galaxy Toolshed and ready to install in any running Galaxy instance.
They are fully workflow compatible and work exactly like any hand-written tool. The user can select input files of the specified type(s) from their
history and edit each of the specified parameters. The tool form will show all the labels and help text supplied when the tool was built. When the tool
is executed, the dependent binary or script will be passed all the i/o files and parameters as specified, and will write outputs to the specified new
history datasets - just like any other Galaxy tool.

## Models for tool execution

The simplest tool model wraps a simple script or Conda dependency package requiring zero parameters, such as filters that take input from STDIN and write to STDOUT.
These can be configured to take STDOUT and write it to a new history item, and to read a user selected history data file on STDIN.

This simple model is illustrated by the Tacrev demonstration tool found in the Galaxy running in the container. It passes a user selected input file on STDIN
to a bash script. The bash script runs the unix tac utility (reverse cat) piped to the unix rev (reverse lines in a text file) utility. It's a one liner:

`tac | rev`

and when run, output from STDOUT is sent to a new history file containing the reversed input text. By reversed, we mean really, truly reversed.

That simple model can be made much more complicated, and can pass inputs and outputs as named or positional parameters,
to allow more complicated scripts or dependent binaries that require:

1. Any number of input data files selected by the user from existing history data
2. Any number of output data files written to the user's history
3. Any number of user supplied parameters. These can be passed as command line arguments to the script or the dependency package. Either
positional or named (argparse) style command line parameter passing can be used.

More complex models can be seen in the Sedtest, Pyrevpos and Pyrevargparse tools illustrating positional and argparse parameter passing.

The most complex demonstration is the Planemo advanced tool tutorial BWA tool. There is one version using a command-override to implement
exactly the same command structure in the Planemo tutorial. A second version uses a bash script and positional parameters to achieve the same
result. Some users may find the bash version more familiar and cleaner but the choice is yours.

## Overview

Steps in building a new Galaxy tool are all conducted through Galaxy running in the docker container:

1. Login to the Galaxy running in the container at http://localhost:8080 using the admin account. They are specified in config/galaxy.yml

2. Start the ToolFactory and fill in the form, providing sample inputs and parameter values to suit the Conda package being wrapped.

3. Execute the tool to create a new XML tool wrapper using the sample inputs and parameter settings for the inbuilt tool test. Planemo is run to generate the outputs
    from the test. The complete toolshed archive is written to the history together with the planemo test report. Optionally the new tool archive can be uploaded
    to the toolshed running in the same container (http://localhost:9009) and then installed inside the Galaxy in the container for further testing.

4. If the test fails, rerun the failed history job and correct errors on the tool form before rerunning until everything works correctly.

![IHello example ToolFactory tool form](files/hello_toolfactory_form.png?raw=true "Part of the Hello world example ToolFactory tool form")


## Planning and building new Galaxy tool wrappers.

It is best to have all the required planning done to wrap any new script or binary before firing up the ToolFactory.
Conda is the only current dependency manager supported. Before starting, at the very least, the user will need
to know the required software package name in Conda and the version to use, how the command line for
the package must be constructed, and there must be sample inputs in the working history for each of the required data inputs
for the package, together with values for every parameter to suit these sample inputs. These are required on the ToolFactory form
for preparing the inbuilt tool test. That test is run using Planemo, as part of the tool generation process.

A new tool is specified by filling in the usual Galaxy tool form.

The form starts with a new tool name. Most tools will need dependency packages and versions
for the executable. Only Conda is currently supported.

If a script is needed, it can be pasted into a text box and the interpreter named. Available system executables
can be used such as bash, or an interpreter such as python, perl or R can be nominated as conda dependencies
to ensure reproducible analyses.

The tool form will be generated from the input data and the user supplied parameters. The command line for the
executable is built using positional or argparse (named e.g. --input_file /foo/baz) style
parameters and is completely dependent on the executable. These can include:

1. Any number of input data sets needed by the executable. Each appears to the tool user on the run form and is included
on the command line for the executable. The tool builder must supply a small representative sample for each one as
an input for the automated tool test.

2. Any number of output data sets generated by the package can be added to the command line and will appear in
the user's history at the end of the job

3. Any number of text or numeric parameters. Each will appear to the tool user on the run form and are included
on the command line to the executable. The tool builder must supply a suitable representative value for each one as
the value to be used for the automated tool test.

Once the form is completed, executing the ToolFactory will build a new XML tool wrapper
including a functional test based on the sample settings and data.

If the Planemo test passes, the tool can be optionally uploaded to the local Galaxy used in the image for more testing.

A local toolshed runs inside the container to allow an automated installation, although any toolshed and any accessible
Galaxy can be specified for this process by editing the default URL and API keys to provide appropriate credentials.

## Generated Tool Dependency management

Conda is used for all dependency management although tools that use system utilities like sed, bash or awk
may be available on job execution nodes. Sed and friends are available as Conda (conda-forge) dependencies if necessary.
Versioned Conda dependencies are always baked-in to the tool and will be used for reproducible calculation.

## Requirements

These are all managed automagically. The ToolFactory relies on galaxyxml to generate tool xml and uses ephemeris and
bioblend to load tools to the toolshed and to Galaxy. Planemo is used for testing and runs in a biocontainer currently at
https://quay.io/fubar2/planemo-biocontainer

This is needed because at present, Planemo seems to have a bug allowing it to leak dependencies back into the calling environment leaving that
environment permanently damaged.  So, it is run completely isolated in a separate container. The docker python SDK is used to manage the
complexities of running docker-in-docker inside the running ToolFactory tool. Trust me - there are complications.

## Caveats

This docker image requires privileged mode so exposes potential security risks if hostile users gain access.
Please, do not run it in any situation where that is a problem - never, ever on a public facing Galaxy server.
On a laptop or workstation should be fine in a non-hostile environment.


## Example generated XML

For the bwa-mem example, a supplied bash script is included as a configfile and so has escaped characters.
```
<tool name="bwatest" id="bwatest" version="0.01">
  <!--Cite: Creating re-usable tools from scripts doi:10.1093/bioinformatics/bts573-->
  <!--Source in git at: https://github.com/fubar2/toolfactory-->
  <!--Created by admin@galaxy.org at 30/11/2020 07:12:10 using the Galaxy Tool Factory.-->
  <description>Planemo advanced tool building sample bwa mem mapper as a ToolFactory demo</description>
  <requirements>
    <requirement version="0.7.15" type="package">bwa</requirement>
    <requirement version="1.3" type="package">samtools</requirement>
  </requirements>
  <configfiles>
    <configfile name="runme"><![CDATA[
REFFILE=\$1
FASTQ=\$2
BAMOUT=\$3
rm -f "refalias"
ln -s "\$REFFILE" "refalias"
bwa index -a is "refalias"
bwa mem -t "2"  -v 1 "refalias" "\$FASTQ"  > tempsam
samtools view -Sb tempsam > temporary_bam_file.bam
samtools sort -o "\$BAMOUT" temporary_bam_file.bam

]]></configfile>
  </configfiles>
  <version_command/>
  <command><![CDATA[bash
$runme
$input1
$input2
$bam_output]]></command>
  <inputs>
    <param optional="false" label="Reference sequence for bwa to map the fastq reads against" help="" format="fasta" multiple="false" type="data" name="input1" argument="input1"/>
    <param optional="false" label="Reads as fastqsanger to align to the reference sequence" help="" format="fastqsanger" multiple="false" type="data" name="input2" argument="input2"/>
  </inputs>
  <outputs>
    <data name="bam_output" format="bam" label="bam_output" hidden="false"/>
  </outputs>
  <tests>
    <test>
      <output name="bam_output" value="bam_output_sample" compare="sim_size" format="bam" delta_frac="0.1"/>
      <param name="input1" value="input1_sample"/>
      <param name="input2" value="input2_sample"/>
    </test>
  </tests>
  <help><![CDATA[

**What it Does**

Planemo advanced tool building sample bwa mem mapper

Reimagined as a bash script for a ToolFactory demonstration


------

Script::

    REFFILE=$1
    FASTQ=$2
    BAMOUT=$3
    rm -f "refalias"
    ln -s "$REFFILE" "refalias"
    bwa index -a is "refalias"
    bwa mem -t "2"  -v 1 "refalias" "$FASTQ"  > tempsam
    samtools view -Sb tempsam > temporary_bam_file.bam
    samtools sort -o "$BAMOUT" temporary_bam_file.bam

]]></help>
</tool>

```






