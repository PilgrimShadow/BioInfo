import com.jgdodson.rosalind.DNAString
import org.scalacheck.Properties
import org.scalacheck.Prop.forAll
import org.scalacheck.Gen

object DNAStringTest extends Properties("DNAString") {

  // Generate non-empty strings from ('A', 'C', 'G', 'T')
  val ACGTGenerator: Gen[String] =
    Gen.containerOf[Vector, Char](Gen.oneOf(DNAString.alphabet)).
      map(_.mkString("")).suchThat(_.length != 0)

  // Generate DNAStrings
  val DNAStrings: Gen[DNAString] = ACGTGenerator.map(DNAString(_))

  property("length") = forAll(ACGTGenerator) { raw =>
    DNAString(raw).length == raw.length
  }

  property("reverse-length") = forAll(DNAStrings) { dna =>
    dna.length == dna.reverse.length
  }

  property("complement-length") = forAll(DNAStrings) { dna =>
    dna.length == dna.complement.length
  }

  property("reverseComplement-length") = forAll(DNAStrings) { dna =>
    dna.length == dna.reverseComplement.length
  }

  property("toRNAString-length") = forAll(DNAStrings) { dna =>
    dna.length == dna.toRNAString.length
  }

  property("complement-gcContent") = forAll(DNAStrings) { dna =>
    dna.GCcontent == dna.complement.GCcontent
  }

  property("complement-atContent") = forAll(DNAStrings) { dna =>
    dna.ATcontent == dna.complement.ATcontent
  }

  property("sum-AT/GC") = forAll(DNAStrings) { dna =>
    val sum = dna.GCcontent + dna.ATcontent
    if (dna.length > 0) sum == 1
    else sum == 0
  }

  property("positive-mass") = forAll(DNAStrings) { dna =>
    dna.mass > 0
  }

  property("positive-masses") = DNAString.masses.forall(_._2 > 0)

  // A, C, G, T
  property("alphabet-size") = DNAString.alphabet.size == 4

}
