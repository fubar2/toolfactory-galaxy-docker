# toolfactory-galaxy-docker

## Docker Galaxy container with the [docker version of the ToolFactory](https://github.com/fubar2/toolfactory_docker) and demonstration tools.

A Galaxy IDE for tool wrappers. Ideal for simple tools. GUI is not ideal for complex ones. Although the ToolFactory should cope, your patience may not.
Turns Galaxy into an integrated development environment for tool builders.
It won't suit every tool or every tool builder but some may find it useful. Perfect for creating simple, custom workflow components requiring system
utilities like sed or awk, or any Conda package and dependencies.

Available from quay.io/fubar2/toolfactory-galaxy-docker but you probably want a
startup script like the ones provided here as samples to set all the ports and volumes.

Updates:<ul><li>Might soon be served in planemo with ``planemo tool_factory --galaxy_root .....``  PR submitted to resurrect the TF</li>
<li>Video *hello world* demonstration https://drive.google.com/file/d/1xpkcVGQ0jRdG78Kt-qLwqeFpE3RnSRsK/view?usp=sharing</li>
</ul>

## Overview

This repository contains the Dockerfile and files needed to build and run the container image available from
quay.io/fubar2/toolfactory-galaxy-docker. It allows the tool builder to build new Galaxy tool wrappers inside Galaxy,
using the specialised ToolFactory tool - a Galaxy tool that generates tool wrappers.

It relies on https://github.com/bgruening/docker-galaxy-stable (20.09 at present) and
the ToolFactory (TF) from https://github.com/fubar2/toolfactory_docker using galaxyxml, ephemeris
and biodocker to do the work of generating and installing new tools. Generated tools are tested with Planemo
from https://github.com/galaxyproject/planemo.

A tool can be built by specifying a Conda or system executable, an optional script and the command line elements needed to run it. The tool will
contain a test and can optionally be installed into the docker Galaxy for further testing.

## Security warning for this container.

*This container must be run in privileged mode, exposing potential security risks*

The TF is a specialised tool that requires a privileged docker container, based on docker-galaxy-stable. Although
the TF tool will only run for administrative tool builders, you are strongly advised not to expose this container on any
public facing production hardware because that potential opportunity for privilege escalation is generally not an acceptable risk.

TF generated tools are just normal Galaxy tool XML wrappers and tests, with no additional security vulnerabilities from the ToolFactory itself.

The TF just makes writing tools easier and uses a familiar Galaxy tool interface. Like any Galaxy tool, TF products *could* be
constructed with malicious code. Outside the privileged docker container used to generate them, generated tools run as normal Galaxy jobs
in an isolated environment where damage will be limited.

## Intended audience and tool builders

Most developers use Planemo on the command line to prepare and test new Galaxy tools. This image
is for tool builders who prefer building and testing new tools *inside the Galaxy UI* !

## Limitations

The TF is flexible enough to generate wrappers for many common scientific packages
but the inbuilt automation will not cope with unusual situations. Tool builders can
supply overrides for two tool XML segments - tests and command and the BWA
example in the supplied samples workflow illustrates their use. It does not deal with
repeated elements or conditional parameters such as allowing a user to choose to see "simple"
or "advanced" parameters (yet) and there will be plenty of packages it just
won't cover - but it's a quick and efficient tool for the other 90% of cases. Perfect for
that bash one liner you need to get that workflow functioning correctly for this
afternoon's demonstration!

## Building the container

Building an image from this Dockerfile is straightforward but most tool builders will probably want to run one of the startup
scripts supplied. They will download the https://quay.io/fubar2/toolfactory_galaxy_docker image or build it, depending
on your preference.

Running the image from a directory with a subdirectory (export/) is recommended for Bjoern's Galaxy docker because
it allows persistence between container restarts. Set this to suit your local needs.

Occasionally, Galaxy fails to launch because the database is not available and after a minute, I just kill it and
start afresh. Usually launches fine then. The script does a thorough cleanup before it starts. Maybe the docker
daemon needs more frequent restarts.

There are two supplied shell scripts with appropriate command lines to either start a fresh clean instance or restart
the image without clearing out the old database and tools:

1. runclean.sh will clear out the export/ directory before running. The first run may take 10-15 minutes to configure and start
Galaxy, a toolshed and install the TF tool with dependencies.

2. rundirty will *not* empty the export/ directory, so can start a previous instance although using docker stop and start is the best way
to preserve a working container.

Both will open ports including 8080/9009. Privileged mode is needed to run the planemo biocontainer as a docker-in-docker image.

Wait until activity, particularly conda and pip, dies down before logging in. A lot happens so expect a 10-15 minute wait on a laptop. When there are no more conda or
other build processes, you should be ok logging in using credentials described at https://github.com/bgruening/docker-galaxy-stable

## Demonstration tools

There are currently a bunch of sample tools built by running a workflow included in the image.
They illustrate the range of models for tool execution that the TF can produce, described in the next section.
Run the workflow, selecting supplied history dataset names as shown to match the filename prompts in each input file field.

![Input file configuration to run the supplied workflow](files/TFWorkflow_setup.png?raw=true "Input file configuration to run the supplied workflow")

The history in which they were all built allows each one to be rerun so you can try re-running them to generating new tools by editing their original build tool form.

See below for generated examples.


## ToolFactory generated tools are ordinary Galaxy tools

A TF generated tool that passes the Planemo test is ready to publish in any Galaxy Toolshed and ready to install in any running Galaxy instance.
They are fully workflow compatible and work exactly like any hand-written tool. The user can select input files of the specified type(s) from their
history and edit each of the specified parameters. The tool form will show all the labels and help text supplied when the tool was built. When the tool
is executed, the dependent binary or script will be passed all the i/o files and parameters as specified, and will write outputs to the specified new
history datasets - just like any other Galaxy tool.

## Models for tool command line construction

The key to turning any software package into a Galaxy tool is the automated construction of a suitable command line.

The TF can build a new tool that will allow the tool user to select input files from their history, set any parameters and when run will send the
new output files to the history as specified when the tool builder completed the form and built the new tool.

That tool can contain instructions to run any Conda dependency or a system executable like bash. Whether a bash script you have written or
a Conda package like bwa, the executable will expect to find settings for input, output and parameters on a command line.

IThese are often passed as "--name value" (argparse style) or in a fixed order (positional style).

The ToolFactory allows either, or for "filter" applications that process input from STDIN and write processed output to STDOUT.

The simplest tool model wraps a simple script or Conda dependency package requiring only input and output files, with no user supplied settings illustrated by
the Tacrev demonstration tool found in the Galaxy running in the ToolFactory docker container. It passes a user selected input file from the current history on STDIN
to a bash script. The bash script runs the unix tac utility (reverse cat) piped to the unix rev (reverse lines in a text file) utility. It's a one liner:

`tac | rev`

The tool building form allows naming zero or more Conda package name(s) and version(s) and the supply of a script to be executed by either a system
executable like ``bash`` or the first of any named Conda dependency package/version. Tacrev uses a tiny bash script shown above and uses the system
bash. Conda bash can be specified if it is important to use the same version consistently for the tool.

On the tool form, the repeat section allowing zero or more input files was set to be a text file to be selected by the tool user and
in the repeat section allowing one or more outputs, a new output file with special value `STDOUT` as the positional parameter, causes the TF to
generate a command to capture STDOUT and send it to the new history file containing the reversed input text.

By reversed, we mean really, truly reversed.

That simple model can be made much more complicated, and can pass inputs and outputs as named or positional parameters,
to allow more complicated scripts or dependent binaries that require:

1. Any number of input data files selected by the user from existing history data
2. Any number of output data files written to the user's history
3. Any number of user supplied parameters. These can be passed as command line arguments to the script or the dependency package. Either
positional or named (argparse) style command line parameter passing can be used.

More complex models can be seen in the Sedtest, Pyrevpos and Pyrevargparse tools illustrating positional and argparse parameter passing.

The most complex demonstration is the Planemo advanced tool tutorial BWA tool. There is one version using a command-override to implement
exactly the same command structure in the Planemo tutorial. A second version uses a bash script and positional parameters to achieve the same
result. Some tool builders may find the bash version more familiar and cleaner but the choice is yours.

## Overview

Steps in building a new Galaxy tool are all conducted through Galaxy running in the docker container:

1. Login to the Galaxy running in the container at http://localhost:8080 using an admin account. They are specified in config/galaxy.yml and
    in the documentation at
    and the ToolFactory will error out and refuse to run for non-administrative tool builders as a minimal protection from opportunistic hostile use.

2. Start the TF and fill in the form, providing sample inputs and parameter values to suit the Conda package being wrapped.

3. Execute the tool to create a new XML tool wrapper using the sample inputs and parameter settings for the inbuilt tool test. Planemo runs twice.
    firstly to generate the test outputs and then to perform a proper test. The completed toolshed archive is written to the history
    together with the planemo test report. Optionally the new tool archive can be uploaded
    to the toolshed running in the same container (http://localhost:9009) and then installed inside the Galaxy in the container for further testing.

4. If the test fails, rerun the failed history job and correct errors on the tool form before rerunning until everything works correctly.

![How it works](files/TFasIDE.png?raw=true "Overview of the ToolFactory and using it as an Integrated Development Environment")

## Planning and building new Galaxy tool wrappers.

It is best to have all the required planning done to wrap any new script or binary before firing up the TF.
Conda is the only current dependency manager supported. Before starting, at the very least, the tool builder will need
to know the required software package name in Conda and the version to use, how the command line for
the package must be constructed, and there must be sample inputs in the working history for each of the required data inputs
for the package, together with values for every parameter to suit these sample inputs. These are required on the TF form
for preparing the inbuilt tool test. That test is run using Planemo, as part of the tool generation process.

A new tool is specified by filling in the usual Galaxy tool form.

The form starts with a new tool name. Most tools will need dependency packages and versions
for the executable. Only Conda is currently supported.

If a script is needed, it can be pasted into a text box and the interpreter named. Available system executables
can be used such as bash, or an interpreter such as python, perl or R can be nominated as conda dependencies
to ensure reproducible analyses.

The hello example uses a bash script to echo a string with a parameter substitution supplied by the user for example

![IHello example ToolFactory tool form](files/hello_toolfactory_form.png?raw=true "Part of the Hello world example ToolFactory tool form")


The tool form will be generated from the input data and the tool builder supplied parameters. The command line for the
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

Once the form is completed, executing the TF will build a new XML tool wrapper
including a functional test based on the sample settings and data.

If the Planemo test passes, the tool can be optionally uploaded to the local Galaxy used in the image for more testing.

A local toolshed runs inside the container to allow an automated installation, although any toolshed and any accessible
Galaxy can be specified for this process by editing the default URL and API keys to provide appropriate credentials.

## Generated Tool Dependency management

Conda is used for all dependency management although tools that use system utilities like sed, bash or awk
may be available on job execution nodes. Sed and friends are available as Conda (conda-forge) dependencies if necessary.
Versioned Conda dependencies are always baked-in to the tool and will be used for reproducible calculation.

## Requirements

These are all managed automagically. The TF relies on galaxyxml to generate tool xml and uses ephemeris and
bioblend to load tools to the toolshed and to Galaxy. Planemo is used for testing and runs in a biocontainer currently at
https://quay.io/fubar2/planemo-biocontainer

## Caveats

This docker image requires privileged mode so exposes potential security risks if hostile tool builders gain access.
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






