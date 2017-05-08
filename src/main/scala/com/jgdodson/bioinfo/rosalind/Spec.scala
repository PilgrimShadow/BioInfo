package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.ProteinString

object Spec {

  def main(args: Array[String]): Unit = {

    if (args.length >= 1) {

      val spectrum = io.Source.fromFile(args(0)).getLines().map(_.toDouble).toVector

      val protein = ProteinString.fromSpectrum(spectrum)

      println(protein.seq)
    }

  }
}
