# Transform4s

[![Build Status](https://travis-ci.org/Sciss/Transform4s.svg?branch=main)](https://travis-ci.org/Sciss/Transform4s)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/de.sciss/transform4s_2.13/badge.svg)](https://maven-badges.herokuapp.com/maven-central/de.sciss/transform4s_2.13)

## statement

This is a library for FFT based transformations. It is a direct translation from Java to Scala of parts
of the [JTransforms](https://github.com/wendykierp/JTransforms) library (C)opyright by Piotr Wendykier. 
All adaptations and Scala translation (C)opyright 2020 by Hanns Holger Ruz. This project is released
under the GNU Lesser General Public License v2.1+.

__Why?__ Because I need a library that can compile both the JVM and JS. In the future, compiling to
Scala Native might also be useful.

__What?__

- Removed all dependencies on other Java libraries (e.g. Apache Commons FastMath, at least for now)
- Removed API unsupported by Scala.js (e.g. Java Logging)
- Removed the Large Array support (not needed in my case,  another Java dependency).

## building

This project builds with sbt against Scala 2.13, 2.12, Dotty (JVM) and 2.13 (Scala.js).

## linking

The following artifact is necessary as dependency:

    libraryDependencies += "de.sciss" %% "transform4s" % v

The current version `v` is `"9.1.0"`
