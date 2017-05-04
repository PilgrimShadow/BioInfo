package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.utils.Utils

object Tree {

  def main(args: Array[String]): Unit = {
    val lines = Utils.readLines(args(0))
    println(lines.head.toInt - lines.length)
  }

}
