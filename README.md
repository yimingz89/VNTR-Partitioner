# VNTR-Partitioner
Partitions VNTRs into individual repeating periods and then aligns them using MUSCLE
## Usage
java -cp [classpath] partitioner.VNTRPartitioner -S [scan file] -R [reference file] -I [VNTR identifier] -O [Root output directory] -U [Updated scan directory]
## Required imports
* SVToolkit.jar
* jaligner.jar
* gatk-utils-3.7.jar
* log4j-1.2.15.jar
* htsjdk-2.19.0-gs.jar
* commons-lang-2.5.jar
* reflections-0.9.9-RC1.jar
