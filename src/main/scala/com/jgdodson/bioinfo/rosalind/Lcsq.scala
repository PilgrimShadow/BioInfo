package com.jgdodson.bioinfo.rosalind

import java.io.{File, PrintWriter}

import com.jgdodson.bioinfo.utils.Utils
import com.jgdodson.bioinfo.DNAString

object Lcsq {

  def main(args: Array[String]): Unit = {

    if (args.length >= 1) {

      val group = Utils.readFastaFile(args(0)).map(fastaEntry => DNAString(fastaEntry._2))

      println()

      // Display the input sequences
      println(s"Sequence 1: (${group(0).length} bp) " + group(0).seq)
      println()
      println(s"Sequence 2: (${group(1).length} bp) " + group(1).seq)
      println()

      val ssm = group(0).longestSharedSplicedMotif(group(1))

      // Display the SSM
      println(s"SSM (${ssm.length} bp): " + ssm)
      println()

      // Display debugging info
      println("Testing: Sequence 1 contains SSM: " + group(0).containsSplicedMotif(DNAString(ssm)))
      println("Testing: Sequence 2 contains SSM: " + group(1).containsSplicedMotif(DNAString(ssm)))

      println()

      // Write the result to a file if second argument is given
      if (args.length >= 2) {
        val pw = new PrintWriter(new File(args(1)))

        pw.write(ssm)

        pw.close()
      }
    }
  }

}
