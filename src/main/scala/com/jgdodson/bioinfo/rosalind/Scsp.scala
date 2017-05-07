package com.jgdodson.bioinfo.rosalind

import java.io.{File, PrintWriter}

import com.jgdodson.bioinfo.DNAString
import com.jgdodson.bioinfo.utils.Utils


object Scsp {

  def main(args: Array[String]): Unit = {
    if (args.length >= 1) {

      val group = Utils.readFastaFile(args(0)).map(fastaEntry => DNAString(fastaEntry._2))

      println()

      // Display the input sequences
      println(s"Sequence 1: (${group(0).length} bp) " + group(0).seq)
      println()
      println(s"Sequence 2: (${group(1).length} bp) " + group(1).seq)
      println()

      val scs = group(0).shortestCommonSuperSequence(group(1))

      // Display the SSM
      println(s"SCS (${scs.length} bp): " + scs)
      println()

      // Display debugging info
      println("Testing: SCS contains Sequence 1: " + DNAString(scs).containsSplicedMotif(group(0)))
      println("Testing: SCS contains Sequence 2: " + DNAString(scs).containsSplicedMotif(group(1)))

      println()

      // Write the result to a file if second argument is given
      if (args.length >= 2) {
        val pw = new PrintWriter(new File(args(1)))

        pw.write(scs)

        pw.close()
      }
    }

  }
}
