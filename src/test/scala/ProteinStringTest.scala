import com.jgdodson.bioinfo.ProteinString
import org.scalacheck.Properties
import org.scalacheck.Prop.forAll
import org.scalacheck.Gen

object ProteinStringTest extends Properties("ProteinStringTest") {

  // Generate non-empty strings from the ProteinString alphabet
  val PeptideStrings: Gen[String] =
    Gen.containerOf[Vector, Char](Gen.oneOf(ProteinString.alphabet)).
      map(_.mkString("")).suchThat(_.length != 0)

  // Generate non-empty strings from the ProteinString alphabet
  val WithoutStopCodon: Gen[ProteinString] =
    Gen.containerOf[Vector, Char](Gen.oneOf(ProteinString.alphabet.filter(_ != 'X'))).
      map(_.mkString("")).suchThat(_.length != 0).map(ProteinString(_))

  // Generate ProteinStrings
  val ProteinStrings: Gen[ProteinString] = PeptideStrings.map(ProteinString(_))

  property("length") = forAll(PeptideStrings) { raw =>
    raw.length == ProteinString(raw).length
  }

  property("reverse-length") = forAll(ProteinStrings) { protein =>
    protein.length == protein.reverse.length
  }

  property("positive-mass") = forAll(WithoutStopCodon) { protein =>
    protein.mass > 0
  }

  property("positive-masses") = ProteinString.masses.forall(_._2 > 0)

  // 20 amino acids + 1 stop codon
  property("alphabet-size") = ProteinString.alphabet.size == 21

}
