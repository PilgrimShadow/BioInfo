package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.ProteinString
import com.jgdodson.bioinfo.utils.Utils


object Edit {

  def main(args: Array[String]): Unit = {

    if (args.length >= 1) {

      val group = Utils.readFastaFile(args(0)).map(fastaEntry => ProteinString(fastaEntry._2))

      // Display the input sequences
      println(s"Sequence 1: (${group(0).length} bp) " + group(0).seq)
      println()
      println(s"Sequence 2: (${group(1).length} bp) " + group(1).seq)
      println()

      println(s"Edit Distance: ${group(0).editDistance(group(1))}")

    }

  }

}
