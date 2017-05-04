package com.jgdodson.bioinfo.rosalind

import java.io.{File, PrintWriter}

import com.jgdodson.bioinfo.utils.Utils
import com.jgdodson.bioinfo.DNAString

object Lcsq {

  def main(args: Array[String]): Unit = {

    if (args.length >= 1) {
      val group = Utils.readFastaFile(args(0)).map(fastaEntry => DNAString(fastaEntry._2))

      val ssm = group(0).longestSharedSplicedMotif(group(1))

      println("Sequence 1: " + group(0).containsSplicedMotif(DNAString(ssm)))
      println("Sequence 2:" + group(1).containsSplicedMotif(DNAString(ssm)))

      if (args.length >= 2) {
        val pw = new PrintWriter(new File(args(1)))

        pw.write(ssm)

        pw.close()
      } else {
        println(ssm)
      }
    }
  }

}
