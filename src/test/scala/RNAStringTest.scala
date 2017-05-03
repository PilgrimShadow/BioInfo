import com.jgdodson.rosalind.{ProteinString, RNAString}
import org.scalacheck.Properties
import org.scalacheck.Prop.forAll
import org.scalacheck.Gen


object RNAStringTest extends Properties("RNAStringTest") {

  // Generate non-empty strings from ('A', 'C', 'G', 'U')
  val ACGUStrings: Gen[String] =
    Gen.containerOf[Vector, Char](Gen.oneOf(RNAString.alphabet)).
      map(_.mkString("")).suchThat(_.length != 0)

  // Generate RNAStrings
  val RNAStrings: Gen[RNAString] = ACGUStrings.map(RNAString(_))

  property("length") = forAll(ACGUStrings) { acgu =>
    acgu.length == RNAString(acgu).length
  }

  property("reverse-length") = forAll(RNAStrings) { rna =>
    rna.length == rna.reverse.length
  }

  property("toDNAString-length") = forAll(RNAStrings) { rna =>
    rna.length == rna.toDNAString.length
  }

  property("valid-translation") = forAll(RNAStrings suchThat (_.length >= 3)) { rna =>
    rna.toProteinString.seq.forall(ch => ProteinString.alphabet.contains(ch))
  }

  property("positive-mass") = forAll(RNAStrings) { rna =>
    rna.mass > 0
  }

  property("positive-masses") = RNAString.masses.forall(_._2 > 0)

  // A,C,G,U
  property("alphabet-size") = RNAString.alphabet.size == 4

}
