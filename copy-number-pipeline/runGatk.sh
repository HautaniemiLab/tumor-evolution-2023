#!/bin/sh

# In principle, JAVA_HOME could be set using PRE scripts. https://anduril.org/userguide/workflows/organizing/#running_scripts_before_and_after_workflow
# However, they also affects Anduril engine and may cause incompatibilities

export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH

/opt/share/gatk-4.1.4.1/gatk "$@"
