lazy val baseName  = "Transform4s"
lazy val baseNameL = baseName.toLowerCase

lazy val projectVersion = "0.1.1"
lazy val mimaVersion    = "0.1.0"

lazy val commonJvmSettings = Seq(
  crossScalaVersions := Seq("3.0.0-M2", "2.13.4", "2.12.12"),
)

// When you want to edit the Java sources in the IDE.
// Disable for publishing!
val ENABLE_JAVA = false

// sonatype plugin requires that these are in global
ThisBuild / version      := projectVersion
ThisBuild / organization := "de.sciss"

lazy val root = crossProject(JSPlatform, JVMPlatform).in(file("."))
  .jvmSettings(commonJvmSettings)
  .settings(
    name               := baseName,
//    version            := projectVersion,
//    organization       := "de.sciss",
    scalaVersion       := "2.13.4",
    description        := "FFT based transforms, such as DFT, DCT, DHT",
    homepage           := Some(url(s"https://github.com/Sciss/${name.value}")),
    licenses           := Seq("BSD 2-Clause" -> url("http://opensource.org/licenses/BSD-2-Clause")),
    mimaPreviousArtifacts := Set("de.sciss" %% baseNameL % mimaVersion),
    unmanagedSourceDirectories in Compile := {
      val all = (unmanagedSourceDirectories in Compile).value
      if (ENABLE_JAVA) all else all.filterNot(_.getPath.contains("java"))
    },
    unmanagedSourceDirectories in Test := {
      val all = (unmanagedSourceDirectories in Test).value
      if (ENABLE_JAVA) all else all.filterNot(_.getPath.contains("java"))
    },
    scalacOptions ++= Seq(
      "-deprecation", "-unchecked", "-feature", "-encoding", "utf8", "-Xlint", "-Xsource:2.13",
    ),
    scalacOptions in (Compile, compile) ++= {
      val jdkGt8  = scala.util.Properties.isJavaAtLeast("9")
      val sv      = scalaVersion.value
      val dot     = isDotty.value // https://github.com/lampepfl/dotty/issues/8634
      val sq0     = if (!dot && jdkGt8) List("-release", "8") else Nil
      if (sv.startsWith("2.12.")) sq0 else "-Wvalue-discard" :: sq0
    }, // JDK >8 breaks API; skip scala-doc
  )
  .jvmSettings(
    libraryDependencies ++= {
      if (!ENABLE_JAVA) Nil else Seq(
        "org.apache.commons"  %   "commons-math3" % "3.5",
        "pl.edu.icm"          %   "JLargeArrays"  % "1.5",
        "junit"               %   "junit"         % "4.13.1" % Test,
      )
    }, 
  )
  .settings(publishSettings)

// ---- publishing ----
lazy val publishSettings = Seq(
  publishMavenStyle := true,
  publishArtifact in Test := false,
  pomIncludeRepository := { _ => false },
  developers := List(
    Developer(
      id    = "p.wendykier",
      name  = "Piotr Wendykier",
      email = "piotr.wendykier@gmail.com",
      url   = url("https://github.com/wendykierp"),
    ),
    Developer(
      id    = "sciss",
      name  = "Hanns Holger Rutz",
      email = "contact@sciss.de",
      url   = url("https://www.sciss.de"),
    )
  ),
  scmInfo := {
    val h = "git.iem.at"
    val a = s"sciss/${name.value}"
    Some(ScmInfo(url(s"https://$h/$a"), s"scm:git@$h:$a.git"))
  },
)

