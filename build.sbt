lazy val baseName  = "Transform4s"
lazy val baseNameL = baseName.toLowerCase

lazy val projectVersion = "0.1.0-SNAPSHOT"
lazy val mimaVersion    = "0.1.0"

lazy val deps = new {
  val test = new {
    val scalaTest = "3.2.3"
  }
}

lazy val commonJvmSettings = Seq(
  crossScalaVersions := Seq("3.0.0-M1", "2.13.3", "2.12.12"),
)

lazy val root = project.in(file(".")) // crossProject(JSPlatform, JVMPlatform).in(file("."))
  .settings(commonJvmSettings)
  .settings(
    name               := baseName,
    version            := projectVersion,
    organization       := "de.sciss",
    scalaVersion       := "2.13.3",
    description        := "FFT based transforms, such as DFT, DCT, DHT",
    homepage           := Some(url(s"https://github.com/Sciss/${name.value}")),
    licenses           := Seq("BSD 2-Clause" -> url("http://opensource.org/licenses/BSD-2-Clause")),
    mimaPreviousArtifacts := Set("de.sciss" %% baseNameL % mimaVersion),
    libraryDependencies ++= Seq(
      "org.apache.commons"  %   "commons-math3" % "3.5",
      "pl.edu.icm"          %   "JLargeArrays"  % "1.5",
      "junit"               %   "junit"         % "4.13.1" % Test,
//      "org.scalatest"       %%  "scalatest"     % deps.test.scalaTest % Test,
    ),
    scalacOptions ++= Seq(
      "-deprecation", "-unchecked", "-feature", "-encoding", "utf8", "-Xlint", "-Xsource:2.13",
    ),
    scalacOptions in (Compile, compile) ++= {
      val jdkGt8  = scala.util.Properties.isJavaAtLeast("9")
      val sv      = scalaVersion.value
      val isDotty = sv.startsWith("3.") // https://github.com/lampepfl/dotty/issues/8634
      val sq0     = (if (!isDotty && jdkGt8) List("-release", "8") else Nil)
      if (sv.startsWith("2.12.")) sq0 else "-Wvalue-discard" :: sq0
    }, // JDK >8 breaks API; skip scala-doc
  )
  .settings(publishSettings)

// ---- publishing ----
lazy val publishSettings = Seq(
  publishMavenStyle := true,
  publishTo := {
    Some(if (isSnapshot.value)
      "Sonatype Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots"
    else
      "Sonatype Releases"  at "https://oss.sonatype.org/service/local/staging/deploy/maven2"
    )
  },
  publishArtifact in Test := false,
  pomIncludeRepository := { _ => false },
  pomExtra := { val n = name.value
<scm>
  <url>git@git.iem.at:sciss/{n}.git</url>
  <connection>scm:git:git@git.iem.at:sciss/{n}.git</connection>
</scm>
<developers>
  <developer>
    <id>p.wendykier</id>
    <name>Piotr Wendykier</name>
    <email>piotr.wendykier@gmail.com</email>
    <timezone>+1</timezone>
  </developer>
  <developer>
    <id>sciss</id>
    <name>Hanns Holger Rutz</name>
    <url>http://www.sciss.de</url>
  </developer>
</developers>
  }
)
