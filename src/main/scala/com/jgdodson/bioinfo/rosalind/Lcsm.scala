package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.utils.Utils

object Lcsm {

  def main(args: Array[String]) = {

    val strings = Utils.readFastaFile(args(0)).map(_._2)
    println(longCommSub(strings.toSet))

  }


  def longCommSub(strings: Set[String]): String = {

    val shortest = strings.minBy(_.length)
    val rest = strings - shortest

    val commonSubs = (0 until shortest.length).foldLeft(Set[String]()) { (acc, i) =>

      val res = ((i + 1) to shortest.length).indexWhere(ind => !rest.forall(_.contains(shortest.substring(i, ind))))

      if (res == 0) acc
      else if (res == -1) acc + shortest.substring(i)
      else acc + shortest.substring(i, i + res)
    }

    commonSubs.maxBy(_.length)
  }
}
