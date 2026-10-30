"""PDB record handling for submission complexes.

The challenge accepts one PDB per compound holding protein + heme + ligand, and the
scorer is strict about things a docking tool has no reason to get right (see
`CLAUDE.md` conventions). This module reads and rewrites those records; it makes no
chemistry decisions.

Why hand-written line handling rather than a parser: every published structure must be
byte-predictable, and the format problems we guard against (records after `END`,
duplicate serials, a heme on the ligand's chain) are properties of the *lines*. A parser
would normalise some of them away silently, which is exactly what hid the heme-dropping
bug in a sibling project.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# Residue names that are neither protein, ligand, nor a cofactor worth keeping.
SOLVENT_RESNAMES = frozenset({"HOH", "WAT", "DOD", "SOL", "TIP3"})

# The ligand residue name the challenge requires.
LIGAND_RESNAME = "LIG"

# Cofactor residue names seen for CYP heme across the PDB. HEM is by far the most
# common; the others appear in some P450 entries and in tool output.
HEME_RESNAMES = frozenset({"HEM", "HEC", "HEB", "HEA", "HDD", "SRM", "DHE"})

STANDARD_AA = frozenset(
    """ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
    MSE SEC PYL ASX GLX UNK HID HIE HIP ASH GLH LYN CYM CYX HSD HSE HSP""".split()
)

_COORD_RECORDS = ("ATOM  ", "HETATM")
_TERMINATORS = ("END", "ENDMDL")


@dataclass
class Atom:
    """One coordinate record, kept as its original line plus parsed fields."""

    line: str

    @property
    def record(self) -> str:
        return self.line[:6]

    @property
    def serial(self) -> str:
        return self.line[6:11]

    @property
    def name(self) -> str:
        return self.line[12:16].strip()

    @property
    def resname(self) -> str:
        return self.line[17:20].strip()

    @property
    def chain(self) -> str:
        return self.line[21:22]

    @property
    def resseq(self) -> str:
        return self.line[22:27].strip()  # resSeq + iCode

    @property
    def element(self) -> str:
        el = self.line[76:78].strip()
        if el:
            return el.upper()
        # Fall back to the atom name's leading alphabetic run. Tools that write short
        # lines omit columns 77-78, and PoseBusters needs the element.
        nm = self.line[12:16].strip()
        return "".join(c for c in nm if c.isalpha())[:1].upper()

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (
            float(self.line[30:38]),
            float(self.line[38:46]),
            float(self.line[46:54]),
        )

    @property
    def residue_key(self) -> tuple[str, str, str]:
        return (self.chain, self.resseq, self.resname)

    def is_protein(self) -> bool:
        """True for standard amino-acid residues only.

        Deliberately ignores the record type: docking tools write ligands as `ATOM`
        often enough that trusting `ATOM` would classify the ligand as protein, and
        the ligand would then be dropped from the submission entirely.
        """
        return self.resname in STANDARD_AA

    def is_solvent(self) -> bool:
        return self.resname in SOLVENT_RESNAMES

    def is_ligand(self) -> bool:
        return self.resname == LIGAND_RESNAME

    def is_heme(self) -> bool:
        return self.resname in HEME_RESNAMES


@dataclass
class Complex:
    """The coordinate records of one complex PDB, split by role.

    Records after a terminator are counted in `after_end`: RDKit's PDB reader stops at
    `END`, so anything past it is invisible to the scorer even though a text search
    would find it.
    """

    path: Path | None = None
    protein: list[Atom] = field(default_factory=list)
    heme: list[Atom] = field(default_factory=list)
    ligand: list[Atom] = field(default_factory=list)
    other: list[Atom] = field(default_factory=list)
    solvent: list[Atom] = field(default_factory=list)
    after_end: list[Atom] = field(default_factory=list)

    @property
    def all_kept(self) -> list[Atom]:
        return self.protein + self.heme + self.ligand + self.other

    def chains(self) -> set[str]:
        return {a.chain for a in self.all_kept}

    def residues(self, atoms: list[Atom] | None = None) -> list[tuple[str, str, str]]:
        seen: list[tuple[str, str, str]] = []
        for a in self.all_kept if atoms is None else atoms:
            if a.residue_key not in seen:
                seen.append(a.residue_key)
        return seen


CIF_SUFFIXES = (".cif", ".mmcif")


def pdb_text(path: str | Path) -> str:
    """Read a structure file as PDB text, converting mmCIF if needed.

    Co-folding tools (Boltz-2, Chai-1) write mmCIF, while the challenge requires PDB.
    gemmi does the conversion; `setup_entities` is needed or the ligand can come back
    without a chain assignment.
    """
    path = Path(path)
    if path.suffix.lower() not in CIF_SUFFIXES:
        return path.read_text(encoding="utf-8", errors="replace")
    import gemmi

    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    return structure.make_pdb_string()


def read_complex(path: str | Path) -> Complex:
    """Parse a complex structure into roles, tracking records stranded after `END`."""
    path = Path(path)
    cx = Complex(path=path)
    ended = False
    for raw in pdb_text(path).splitlines():
        rec = raw[:6]
        if raw[:6].rstrip() in _TERMINATORS or raw[:3] == "END":
            ended = True
            continue
        if rec not in _COORD_RECORDS or len(raw) < 54:
            continue
        atom = Atom(raw.ljust(80))
        # Classified by residue name, in priority order. Anything unrecognised goes to
        # `other` rather than being guessed at, so the caller sees it (a docking tool's
        # ligand usually lands there, under MOL/UNL/<code>).
        if ended:
            cx.after_end.append(atom)
        elif atom.is_ligand():
            cx.ligand.append(atom)
        elif atom.is_heme():
            cx.heme.append(atom)
        elif atom.is_solvent():
            cx.solvent.append(atom)
        elif atom.is_protein():
            cx.protein.append(atom)
        else:
            cx.other.append(atom)
    return cx


def _set(line: str, start: int, value: str, width: int) -> str:
    """Overwrite a fixed-width PDB column span, right-justified."""
    return line[:start] + value.rjust(width)[-width:] + line[start + width :]


def unique_atom_names(atoms: list[Atom]) -> list[Atom]:
    """Make atom names unique within each residue, changing nothing that already is.

    The validator rejects duplicate atom names in one residue because OST's
    symmetry-corrected matching needs them to tell atoms apart, and docking output
    routinely writes every carbon as `C`. Where a rename is needed, names become
    element + running index per residue (C1, C2, N1 ...).

    Residues whose names are already unique are left **untouched**, which matters far
    more than it looks: a cofactor carries standard CCD atom names (HEM's `FE`, `NA`,
    `CHA`, `C1A` ...) and the scorer matches the model's cofactor to the reference's by
    those names. Renaming a correctly named heme to `FE1`/`N1`/`C1` makes every one of
    its atoms unmatchable, which silently deletes all ligand-heme contacts from LDDT-PLI
    and drops the heme out of the binding-site LDDT-LP. Renaming unconditionally is
    therefore not a harmless normalisation; it is a scoring bug.
    """
    by_residue: dict[tuple[str, str, str], list[Atom]] = {}
    for a in atoms:
        by_residue.setdefault(a.residue_key, []).append(a)
    needs_rename = set()
    for key, residue_atoms in by_residue.items():
        names = [a.name for a in residue_atoms]
        if len(names) != len(set(names)) or any(not n for n in names):
            needs_rename.add(key)

    out: list[Atom] = []
    counters: dict[tuple[tuple[str, str, str], str], int] = {}
    for a in atoms:
        if a.residue_key not in needs_rename:
            out.append(a)
            continue
        el = a.element or "X"
        key = (a.residue_key, el)
        counters[key] = counters.get(key, 0) + 1
        name = f"{el}{counters[key]}"
        # Columns 13-16: a one-letter element starts in column 14 by convention.
        field_ = f" {name:<3s}" if len(el) == 1 and len(name) <= 3 else f"{name:<4s}"
        out.append(Atom(_set(a.line, 12, field_[:4].ljust(4), 4)))
    return out


def renumber(atoms: list[Atom], start: int = 1) -> list[Atom]:
    """Give every record a unique serial, in order, starting at `start`."""
    return [
        Atom(_set(a.line, 6, str(start + i), 5)) for i, a in enumerate(atoms)
    ]


def set_chain_and_resseq(atoms: list[Atom], chain: str, resseq: int) -> list[Atom]:
    """Put every record on one chain and residue number (for heme or ligand)."""
    out = []
    for a in atoms:
        line = _set(a.line, 21, chain, 1)
        line = _set(line, 22, str(resseq), 4)
        line = line[:26] + " " + line[27:]  # clear the insertion code
        out.append(Atom(line))
    return out


# Atom names that denote the iron inside a heme, by CCD convention.
HEME_IRON_NAMES = frozenset({"FE", "FE1", "FE2"})


def fix_heme_iron_element(atoms: list[Atom]) -> tuple[list[Atom], int]:
    """Correct an iron whose element column was derived from its name's first letter.

    Some exports write the element by taking the atom name's leading character, which is
    right for `NA` -> N and `CAA` -> C but turns the heme iron's `FE` into **F**. The
    submitted file would then declare a fluorine where the iron belongs, and every
    consumer that reads the element — the scorer, PoseBusters' covalent-radius check,
    RDKit — would believe it.

    Only an atom named FE inside a heme residue is touched, and only when its element is
    not already FE, so this cannot reinterpret a genuine fluorine.

    Returns the atoms and how many were corrected.
    """
    out: list[Atom] = []
    fixed = 0
    for a in atoms:
        if a.is_heme() and a.name.upper() in HEME_IRON_NAMES and a.element != "FE":
            out.append(Atom(_set(a.line, 76, "FE", 2)))
            fixed += 1
        else:
            out.append(a)
    return out, fixed


def heme_iron(cx: "Complex") -> Atom | None:
    """The heme iron, found by element or, failing that, by its CCD atom name.

    The name fallback exists because of the mis-declared-element case above: a run whose
    iron is labelled F must still be measurable, or every coordination metric silently
    reports nothing.
    """
    for atom in cx.heme:
        if atom.element == "FE":
            return atom
    for atom in cx.heme:
        if atom.name.upper() in HEME_IRON_NAMES:
            return atom
    return None


def set_resname(atoms: list[Atom], resname: str, as_hetatm: bool = True) -> list[Atom]:
    """Rename a residue and, by default, force its records to HETATM.

    The challenge requires the ligand residue to be named exactly `LIG`, and a
    non-polymer residue belongs in HETATM records. A docking tool that emits the
    ligand as `MOL`/`UNL`, or as `ATOM`, fails the validator on both counts.
    """
    out = []
    for a in atoms:
        line = _set(a.line, 17, resname.ljust(3)[:3], 3)
        if as_hetatm:
            line = "HETATM" + line[6:]
        out.append(Atom(line))
    return out


def atom_line(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    xyz: tuple[float, float, float],
    element: str,
    record: str = "HETATM",
) -> Atom:
    """Build one coordinate record at the correct PDB columns.

    Hand-formatting these by f-string is a reliable way to end up one column out, which
    does not fail loudly: the coordinates still parse while the residue name and element
    land in the wrong fields, and the element is what PoseBusters needs to perceive the
    molecule. Built once here, with a test, rather than at each call site.

    Columns, 1-indexed per the PDB spec: 1-6 record, 7-11 serial, 13-16 name,
    18-20 resName, 22 chainID, 23-26 resSeq, 31-38/39-46/47-54 x/y/z, 55-60 occupancy,
    61-66 tempFactor, 77-78 element.
    """
    # A one- or two-character element is conventionally right-justified in columns 13-14
    # of the four-column name field, which keeps e.g. CL distinct from a carbon named CL.
    label = name.strip()[:4]
    padded = f" {label:<3}" if len(label) < 4 else label
    return Atom(
        f"{record:<6}"
        f"{serial:>5d} "
        f"{padded}"
        f" {resname.strip()[:3]:>3} "
        f"{chain.strip()[:1] or ' '}"
        f"{resseq:>4d}    "
        f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}"
        f"{1.0:>6.2f}{0.0:>6.2f}"
        f"{'':>10}"
        f"{element.strip().upper()[:2]:>2}"
    )


def write_complex(
    path: str | Path,
    protein: list[Atom],
    heme: list[Atom],
    ligand: list[Atom],
    header: list[str] | None = None,
) -> Path:
    """Write a submission-ready complex: protein, then heme, then ligand, then `END`.

    The ligand is renamed `LIG` and the heme `HEM`, both as HETATM. Chains are fixed as
    A (protein), B (heme), C (ligand) and serials run 1..N across the file, so none of
    the validator's numbering checks can trip. No CONECT records are written: the scorer
    derives ligand bonds from the SMILES template, and a stale CONECT block that points
    at renumbered serials is worse than none.
    """
    path = Path(path)
    if heme:
        heme = set_resname(unique_atom_names(heme), "HEM")
        heme = set_chain_and_resseq(heme, "B", 900)
    ligand = set_resname(unique_atom_names(ligand), LIGAND_RESNAME)
    ligand = set_chain_and_resseq(ligand, "C", 901)
    body = renumber(protein + heme + ligand)
    lines = list(header or [])
    lines += [a.line.rstrip() for a in body]
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    # PDB records are ASCII, but be explicit so output never depends on the locale.
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
