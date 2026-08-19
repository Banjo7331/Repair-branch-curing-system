"""The vocabulary of a mechanism call.

Nothing in this module decides anything. It defines what a variant, a
transcript and a gene *are* for the purposes of asking why a variant causes
disease, and what an answer to that question looks like.

The rules that turn the first into the second do not live in Python at all —
they live in ``rules/mechanism-*.yaml``, because they are clinical claims from
the literature and the person best placed to check them is a geneticist, not a
programmer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class RepairbenchError(Exception):
    """Base class for domain rule violations."""


class InvalidTranscriptError(RepairbenchError):
    """The transcript structure cannot support a prediction."""


class InconsistentGeneError(RepairbenchError):
    """Curated gene context contradicts itself."""


# --------------------------------------------------------------------------
# What a variant does to the sequence
# --------------------------------------------------------------------------


class Consequence(StrEnum):
    """Molecular consequence, in the Sequence Ontology's vocabulary.

    Deliberately coarse. Finer distinctions exist, but every one of them that
    this module would act on is already captured by the transcript arithmetic
    or by the splice prediction, and a longer list would only invite rules that
    look precise while resting on an annotator's judgement call.
    """

    NONSENSE = "stop_gained"
    FRAMESHIFT = "frameshift_variant"
    SPLICE_ACCEPTOR = "splice_acceptor_variant"
    SPLICE_DONOR = "splice_donor_variant"
    SPLICE_REGION = "splice_region_variant"
    MISSENSE = "missense_variant"
    INFRAME_DELETION = "inframe_deletion"
    INFRAME_INSERTION = "inframe_insertion"
    START_LOST = "start_lost"
    STOP_LOST = "stop_lost"
    SYNONYMOUS = "synonymous_variant"
    INTRONIC = "intron_variant"

    @property
    def is_predicted_null(self) -> bool:
        """Does this consequence introduce a premature termination codon?

        Splice variants are excluded even though they often destroy the
        protein: whether they do depends on what the mis-splicing produces, and
        that is a separate question with its own rules.
        """
        return self in {Consequence.NONSENSE, Consequence.FRAMESHIFT, Consequence.START_LOST}

    @property
    def is_splice_affecting(self) -> bool:
        return self in {
            Consequence.SPLICE_ACCEPTOR,
            Consequence.SPLICE_DONOR,
            Consequence.SPLICE_REGION,
        }

    @property
    def preserves_reading_frame(self) -> bool:
        return self in {
            Consequence.MISSENSE,
            Consequence.INFRAME_DELETION,
            Consequence.INFRAME_INSERTION,
            Consequence.SYNONYMOUS,
        }


class Zygosity(StrEnum):
    """How many working copies the patient actually has.

    Absent from the first version of this package, and its absence produced a
    wrong answer rather than a missing one: every modality that works by raising
    output from the intact allele was offered for a hemizygous boy with Duchenne
    muscular dystrophy, who has no intact allele to raise.

    That is the shape of the mistake this enum exists to prevent. A mechanism is
    a property of a variant and a gene; how much functional product the patient
    has left is a property of the *patient*, and half the interventions in M6
    depend on the second rather than the first.
    """

    HETEROZYGOUS = "heterozygous"
    HOMOZYGOUS = "homozygous"
    HEMIZYGOUS = "hemizygous"
    COMPOUND_HETEROZYGOUS = "compound_heterozygous"
    UNKNOWN = "unknown"

    @property
    def leaves_a_wild_type_allele(self) -> bool | None:
        """Is there an unaffected copy left to work with?

        ``None`` for unknown, and the distinction is load-bearing: a rule that
        treated unknown as "no" would rule out real options on missing data,
        and one that treated it as "yes" would offer options that may not exist.
        Unknown raises a caveat instead of deciding either way.
        """
        if self is Zygosity.HETEROZYGOUS:
            return True
        if self is Zygosity.UNKNOWN:
            return None
        # Homozygous, hemizygous and compound heterozygous all mean the same
        # thing here: no copy of this gene is free of pathogenic variation.
        return False


# --------------------------------------------------------------------------
# Gene-level context
# --------------------------------------------------------------------------


class DosageScore(StrEnum):
    """ClinGen dosage sensitivity ratings, kept in their own vocabulary.

    Three of these are not points on a scale, and flattening any of them into
    the others loses a distinction the resolver needs.

    ``UNLIKELY`` is a *finding*: somebody looked and concluded that dosage
    sensitivity does not explain this gene. ``NO_EVIDENCE`` is a *gap in the
    evidence*: ClinGen curated the gene and found nothing to support it yet. And
    ``NOT_EVALUATED`` is a gap in the *curation*: nobody has looked at all.

    The third was added after pointing this package at the real ClinGen list.
    *KRT14* and *SMN1* are simply not on it — and the default here was
    ``NO_EVIDENCE``, so both came out carrying a claim ClinGen has never made.
    A default that silently attributes a finding to a source is the exact
    failure this project's provenance layer exists to prevent, and it was
    sitting in the model the whole time.
    """

    NOT_EVALUATED = "not_evaluated"
    NO_EVIDENCE = "no_evidence"
    LITTLE_EVIDENCE = "little_evidence"
    EMERGING_EVIDENCE = "emerging_evidence"
    SUFFICIENT_EVIDENCE = "sufficient_evidence"
    AUTOSOMAL_RECESSIVE = "autosomal_recessive"
    UNLIKELY = "dosage_sensitivity_unlikely"

    @property
    def supports_haploinsufficiency(self) -> bool:
        return self in {DosageScore.EMERGING_EVIDENCE, DosageScore.SUFFICIENT_EVIDENCE}

    @property
    def refutes_haploinsufficiency(self) -> bool:
        return self in {DosageScore.UNLIKELY, DosageScore.AUTOSOMAL_RECESSIVE}

    @property
    def is_curated(self) -> bool:
        """Did anybody look?

        Every other property here answers "what did they find". This one
        answers whether the question was asked, and the two are only the same
        if you assume a curation exists — which for most of the genome it does
        not: ClinGen has evaluated a few thousand genes, not twenty thousand.
        """
        return self is not DosageScore.NOT_EVALUATED


class Imprinting(StrEnum):
    """Whether an intact copy of the gene is present but switched off.

    Modelled as a property of the gene rather than as a mechanism, on purpose.
    Angelman syndrome is loss of function of *UBE3A* — the imprinting does not
    change what went wrong, it changes what can be done about it, because an
    intact silenced copy is a therapeutic target that most loss-of-function
    genes do not offer.

    ``X_INACTIVATED`` is not imprinting and is included anyway. In a
    heterozygous female, X inactivation leaves the wild-type allele intact but
    silent in about half of cells — biologically a different mechanism,
    therapeutically the same opportunity, and a real research route in Rett
    syndrome. Keeping it out of this enum on taxonomic grounds would have cost
    the model the one modality that case actually has.
    """

    NONE = "not_imprinted"
    MATERNALLY_EXPRESSED = "maternally_expressed"
    PATERNALLY_EXPRESSED = "paternally_expressed"
    X_INACTIVATED = "x_inactivated"

    @property
    def is_imprinted(self) -> bool:
        """Retained name; reads as "an intact copy is silenced"."""
        return self is not Imprinting.NONE


@dataclass(frozen=True, slots=True)
class MissenseDistribution:
    """Where pathogenic variation sits in a gene.

    An indirect signal, but a strong one: a gene whose pathogenic missense
    variants cluster in one domain while its truncating variants are absent is
    a gene where disease comes from the protein doing something new, not from
    there being less of it.
    """

    pathogenic_missense_total: int = 0
    pathogenic_missense_in_hotspot: int = 0
    pathogenic_truncating_total: int = 0

    def __post_init__(self) -> None:
        if self.pathogenic_missense_in_hotspot > self.pathogenic_missense_total:
            raise InconsistentGeneError(
                "more hotspot missense variants than missense variants"
            )

    @property
    def clustering_ratio(self) -> float:
        if self.pathogenic_missense_total == 0:
            return 0.0
        return self.pathogenic_missense_in_hotspot / self.pathogenic_missense_total

    @property
    def counted(self) -> int:
        """How many pathogenic variants went into this at all.

        Nothing reads a ratio without first reading this. A ratio computed from
        six reports and one computed from six hundred are the same number and
        different evidence, and the rule file has a minimum for exactly that
        reason.
        """
        return self.pathogenic_missense_total + self.pathogenic_truncating_total

    @property
    def truncating_fraction(self) -> float | None:
        """What share of pathogenic variation in this gene destroys the product.

        This replaced a boolean — *are there any truncating variants* — which
        real data killed on its first contact with it. *PIK3CA* is the gene the
        rule was written for, and ClinVar reports four truncating variants in it
        against seventy-six missense: enough to make "any" true and nowhere near
        enough to mean truncation is a mechanism there. A curated database of a
        hundred submissions will always contain a few rows that disagree with
        the rest, and a boolean hands each of them a veto.

        ``None`` when nothing was counted, which is not the same as zero. A gene
        nobody has submitted variants for must not read as a gene where
        truncation causes no disease — that is half of the gain-of-function
        argument, offered for free.
        """
        if self.counted == 0:
            return None
        return self.pathogenic_truncating_total / self.counted

    def __str__(self) -> str:
        """The three counts and the ratio, in a line a report can print.

        The ratio is spelled out because it is what the clustering rule reads,
        and a reader checking whether a gain-of-function call was reasonable
        should not have to do the division to see how near a threshold it fell.
        An empty distribution says so rather than printing zeroes, which read
        as a measurement that came back negative.
        """
        if self.pathogenic_missense_total == 0 and self.pathogenic_truncating_total == 0:
            return "nothing counted"
        return (
            f"{self.pathogenic_missense_in_hotspot}/{self.pathogenic_missense_total} "
            f"missense in the densest window ({self.clustering_ratio:.0%}), "
            f"{self.pathogenic_truncating_total} truncating"
        )


@dataclass(frozen=True, slots=True)
class Gene:
    """Curated and population-derived context for one gene."""

    symbol: str
    #: What ClinGen says about losing a copy, and about gaining one. The
    #: default is *not evaluated* rather than *no evidence*: most of the genome
    #: has never been curated, and defaulting to a score would put a claim in
    #: ClinGen's mouth for every gene nobody has looked at.
    haploinsufficiency: DosageScore = DosageScore.NOT_EVALUATED
    triplosensitivity: DosageScore = DosageScore.NOT_EVALUATED
    #: gnomAD loss-of-function observed/expected upper bound. Low means the
    #: population tolerates loss of this gene poorly.
    loeuf: float | None = None
    #: Does the product assemble with copies of itself or with partners? A
    #: mutant subunit can only poison a complex it is part of, so this gates
    #: every dominant-negative rule.
    forms_multimer: bool = False
    #: Published observation that null carriers are milder than missense
    #: carriers. The classic dominant-negative signature: if having none is
    #: better than having a broken one, the broken one is doing harm.
    truncating_variants_are_milder: bool = False
    imprinting: Imprinting = Imprinting.NONE
    #: Whether the silenced parental copy is free of pathogenic variants —
    #: without which reactivation has nothing to reactivate.
    silenced_allele_intact: bool = False
    distribution: MissenseDistribution = field(default_factory=MissenseDistribution)
    #: An expert determination, when one exists. It outranks anything inferred.
    curated_mechanism: str | None = None
    curated_mechanism_source: str | None = None

    def __post_init__(self) -> None:
        if not self.symbol:
            raise InconsistentGeneError("gene symbol is empty")
        if self.curated_mechanism and not self.curated_mechanism_source:
            raise InconsistentGeneError(
                f"{self.symbol}: a curated mechanism without a source cannot outrank inference"
            )
        if self.silenced_allele_intact and not self.imprinting.is_imprinted:
            raise InconsistentGeneError(
                f"{self.symbol}: an intact silenced allele was claimed for a non-imprinted gene"
            )


# --------------------------------------------------------------------------
# The answer
# --------------------------------------------------------------------------


class Mechanism(StrEnum):
    """How a variant causes disease — not what it does to the sequence.

    This distinction is why the module exists. Two nonsense variants in one
    gene share a consequence and can have opposite mechanisms; a mechanism,
    unlike a consequence, tells you what a therapy would have to accomplish.
    """

    LOSS_OF_FUNCTION = "loss_of_function"
    """Less product than the cell needs. The question is how to supply more."""

    GAIN_OF_FUNCTION = "gain_of_function"
    """The product does something it should not. More normal protein does not help."""

    DOMINANT_NEGATIVE = "dominant_negative"
    """The mutant interferes with the normal one — the mechanism where adding a
    working copy can make things worse."""

    SPLICING_DISRUPTION = "splicing_disruption"
    """The transcript is assembled wrongly; the protein-level consequence follows."""

    UNDETERMINED = "undetermined"
    """The evidence does not settle it. A first-class answer, not a failure."""

    @property
    def is_determined(self) -> bool:
        return self is not Mechanism.UNDETERMINED

    @property
    def tolerates_gene_addition(self) -> bool:
        """Is supplying a working copy mechanistically coherent?

        Before any question of vector, tissue or dose — and false for the two
        mechanisms where more normal protein does not address the problem, one
        of which it can actively worsen.
        """
        return self is Mechanism.LOSS_OF_FUNCTION


class Confidence(StrEnum):
    """How firmly the mechanism was established."""

    ESTABLISHED = "established"
    """An expert curation said so."""

    PROBABLE = "probable"
    """Independent rules agree and nothing of weight contradicts them."""

    POSSIBLE = "possible"
    """Something fired, but not enough to lean on."""

    NONE = "none"
    """The rules did not converge."""


class Strength(StrEnum):
    """The weight one rule carries."""

    DECISIVE = "decisive"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"

    @property
    def points(self) -> int:
        return {
            Strength.DECISIVE: 100,
            Strength.STRONG: 3,
            Strength.MODERATE: 2,
            Strength.SUPPORTING: 1,
        }[self]


@dataclass(frozen=True, slots=True)
class Evidence:
    """One rule that fired, what it argues for, and why."""

    rule_id: str
    supports: Mechanism
    strength: Strength
    because: str
    citation: str = ""

    @property
    def label(self) -> str:
        return f"{self.rule_id}({self.strength.value[:3]})"


@dataclass(frozen=True, slots=True)
class Feasibility:
    """What the mechanism and the transcript together make possible.

    These are necessary conditions, never sufficient ones. A true flag means
    "not ruled out on mechanistic grounds", which is a far weaker claim than
    "this would work" — the field names are chosen so a reader is not tempted
    to hear the stronger one.
    """

    gene_addition_coherent: bool = False
    fits_viral_payload: bool = False
    exon_skipping_preserves_frame: bool | None = None
    silenced_allele_available: bool = False
    allele_specific_silencing_indicated: bool = False
    notes: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class MechanismCall:
    """The module's output: a mechanism, its confidence, the rules that produced
    it, and what that mechanism makes possible.

    There is no way to construct one without its evidence, and that is
    deliberate — the same discipline an ACMG classifier keeps about criteria.
    """

    gene: str
    transcript: str
    mechanism: Mechanism
    confidence: Confidence
    evidence: tuple[Evidence, ...]
    feasibility: Feasibility
    #: Rules that argued for a different mechanism. A call with conflicts is
    #: still a call; hiding them would be the problem.
    conflicts: tuple[Evidence, ...] = ()
    #: Digest of the rule file that produced this call. Same discipline as a
    #: knowledge snapshot pin: a verdict that cannot name the rules behind it
    #: cannot be reproduced or compared with a later one.
    ruleset_version: str = ""

    @property
    def needs_review(self) -> bool:
        return (
            not self.mechanism.is_determined
            or self.confidence is Confidence.NONE
            or bool(self.conflicts)
        )

    def summary(self) -> str:
        line = f"{self.gene} ({self.transcript})  {self.mechanism}  [{self.confidence}]"
        if self.evidence:
            line += "  " + ", ".join(e.label for e in self.evidence)
        if self.conflicts:
            line += f"  ⚠ {len(self.conflicts)} conflicting"
        return line
