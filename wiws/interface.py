import os
import logging
import shutil
from tempfile import mkdtemp
from enum import Enum
from typing import List, Tuple, Optional, Dict
import subprocess
from uuid import uuid4

import flask


_log = logging.getLogger(__name__)


bp = flask.Blueprint('api', __name__, url_prefix='/api')


# Description of WHAT IF web services


def run_whatif(command: str, pdbid: str, use_prodrug: Optional[bool] = False):

    from wiws import settings as config

    script = "INIALL\n"
    script += "SETWIF 368 1\n"
    if use_prodrug:
        script += "SETWIF 1052 1\n"
    else:
        script += "SETWIF 1052 0\n"

    if len(pdbid) == 4:
        pdb_path = os.path.join(config.PDB_DIRECTORY_PATH, f"pdb{pdbid.lower()}.ent")
    else:
        pdb_path = os.path.join(config.PDB_UPLOAD_DIRECTORY_PATH, f"{pdbid.lower()}.pdb")

    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(pdb_path)

    script += f"GETMOL {pdb_path} {pdbid}\n"
    script += command + '\n'
    script += "END y\n"

    tmp_dir = mkdtemp()
    try:
        completed_process = subprocess.run([config.WHATIF_EXECUTABLE_PATH],
                                            cwd=tmp_dir,
                                            input=script.encode('ascii'),
                                            capture_output=True, timeout=60.0)
    finally:
        shutil.rmtree(tmp_dir)

    if completed_process.returncode != 0:
        raise RuntimeError(f"whatif returned status {completed_process.returncode}:\n {completed_process.stderr.decode('ascii')}")

    result = completed_process.stdout.decode('ascii').split('\n')
    start_index = None
    end_index = None
    for index, line in enumerate(result):
        if line == '****ENDHEADER****':
            start_index = index
        if line == '****ENDOUTPUT****':
            end_index = index

    result = result[start_index + 1: end_index]

    _log.debug(result)

    return result


class Residue:
    def __init__(self,
                 number: int,
                 chain: str,
                 type_: str,
                 pdb_number: int,
                 insertion_code: str,
                 model_number: Optional[int] = None):

        self.number = number
        self.chain = chain
        self.type = type_
        self.pdb_number = pdb_number
        self.insertion_code = insertion_code
        self.model_number = model_number


class Atom:
    def __init__(self,
                 residue: Residue,
                 atom_type: str,
                 alternate_atom: str):

        self.residue = residue
        self.atom_type = atom_type
        self.alternate_atom = alternate_atom


class FullAtom:
    def __init__(self,
                 PDBx_atom_site: int,
                 PDBx_Cartn_x: float,
                 PDBx_Cartn_y: float,
                 PDBx_Cartn_z: float,
                 PDBx_B_iso_or_equiv: float,
                 PDBx_auth_asym_id: str,
                 PDBx_auth_atom_id: str,
                 PDBx_occupancy: float,
                 PDBx_type_symbol: str):

        self.PDBx_atom_site = PDBx_atom_site
        self.PDBx_Cartn_x = PDBx_Cartn_x
        self.PDBx_Cartn_y = PDBx_Cartn_y
        self.PDBx_Cartn_z = PDBx_Cartn_z
        self.PDBx_B_iso_or_equiv = PDBx_B_iso_or_equiv
        self.PDBx_auth_asym_id = PDBx_auth_asym_id
        self.PDBx_auth_atom_id = PDBx_auth_atom_id
        self.PDBx_occupancy = PDBx_occupancy
        self.PDBx_type_symbol = PDBx_type_symbol


class FullResidue:
    def __init__(self, number: int, type_: str, pdb_number: int, insertion_code: str, chain: str, model_number: int, atoms: List[FullAtom]):
        self.number = number
        self.type = type_
        self.pdb_number = pdb_number
        self.insertinon_code = insertion_code
        self.chain = chain
        self.model_number = model_number
        self.atoms = atoms


class Coordinate:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class Hue(Enum):
    BLUE = 0  # 360
    PURPLE = 60
    RED = 120
    ORANGE = 150
    YELLOW = 180
    GREEN = 240
    LIME = 30


class Dot:
    def __init__(self, coordinate: Coordinate, hue: int):
        self.coordinate = coordinate
        self.hue = hue

class Line:
    def __init__(self, coordinate_from: Coordinate, coordinate_to: Coordinate, hue: int):
        self.coordinate_from = coordinate_from
        self.coordinate_to = coordinate_to
        self.hue = hue

class SurfaceDot:
    def __init__(self, atom: Atom, dots: List[Dot]):
        self.atom = atom
        self.dots = dots


class Rotamer:
    def __init__(self, lines: List[Line], hue: int):
        self.lines = lines
        self.hue = hue


class TorsAng:
    def __init__(self, phi: float, psi: float, omega: float, chi: List[float]):
        self.phi = phi
        self.psi = psi
        self.omega = omega
        self.chi = chi


class TorsAngBB:
    def __init__(self, phi: float, psi: float, omega: float):
        self.phi = phi
        self.psi = psi
        self.omega = omega


class ResidueTorsionInfo:
    def __init__(self, residue: Residue, angles: TorsAng):
        self.residue = residue
        self.angles = angles


class ResidueTorsionInfoBackbone:
    def __init__(self, residue: Residue, angles: TorsAngBB):
        self.residue = residue
        self.angles = angles


class ResidueAngleInfoTau:
    def __init__(self, residue: Residue, tau: float):
        self.residue = residue
        self.tau = tau


class SymmetryContactInfo:
    def __init__(self, residue: Residue, contact_count: int):
        self.residue = residue
        self.contact_count = contact_count


class AtomAccessibilityInfo:
    def __init__(self, atom: Atom, accessibility: float):
        self.atom = atom
        self.accessibility = accessibility


class AtomAccessibilityInfoPlus:
    def __init__(self, atom: Atom, accessibility: float, status: bool):
        self.atom = atom
        self.accessibility = accessibility
        self.status = status


class ResidueAccessibilityInfo:
    def __init__(self, residue: Residue, access_bb: float, access_sch: float, access_tot: float, status: bool):
        self.residue = residue
        self.access_bb = access_bb
        self.access_sch = access_sch
        self.access_tot = access_tot
        self.status = status


class EntityAccessibilityInfo:
    def __init__(self, residue: Residue, access: float, status: bool):
        self.residue = residue
        self.access = access
        self.status = status


class VacuumAccessibilityInfo:
    def __init__(self, residue: Residue, normal_acc: float, vacuum_acc: float):
        self.residue = residue
        self.normal_acc = normal_acc
        self.vacuum_acc = vacuum_acc


class CysteineBridgeInfo:
    def __init__(self, cysteine1: Residue, cysteine2: Residue):
        self.cysteine1 = cysteine1
        self.cysteine2 = cysteine2


class CysteineFreeInfo:
    def __init__(self, cysteine: Residue):
        self.cysteine = cysteine


class CysteineMetalInfo:
    def __init__(self, cysteine: Residue, metal: Residue):
        self.cysteine = cysteine
        self.metal = metal


class CysteineTorsionsInfo:
    def __init__(self, cysteines: Tuple[Residue, Residue],
                 angles: Tuple[float, float, float, float, float]):
        self.cysteines = cysteines
        self.angles = angles


class HydrogenBondInfo:
    def __init__(self, atoms: Tuple[Atom, Atom],
                 hydrogenbond_parameters: Tuple[float, float, float, float]):
        self.atoms = atoms
        self.hydrogenbond_parameters = hydrogenbond_parameters


class BumpInfo:
    def __init__(self, atoms: Tuple[Atom, Atom], bump_severity: float):
        self.atoms = atoms
        self.bump_severity = bump_severity


class LigandContactInfo:
    def __init__(self, atoms: Tuple[Atom, Atom], contact_distance: float):
        self.atoms = atoms
        self.contact_distance = contact_distance


class NucleicAcidContactInfo:
    def __init__(self, atoms: Tuple[Atom, Atom], contact_distance: float):
        self.atoms = atoms
        self.contact_distance = contact_distance


class SaltBridgeInfo:
    def __init__(self, atoms: Tuple[Atom, Atom], distance: float):
        self.atoms = atoms
        self.distance = distance


class LikelyRotamerInfo:
    def __init__(self, residue: Residue, rotamer_table: List[float], quality_estimate: float):
        self.residue = residue
        self.rotamer_table = rotamer_table
        self.quality_estimate = quality_estimate


class AccessAndSymmetryInfo:
    def __init__(self, residue: Residue, acc_sym_values: Tuple[float, float, float, float, float, float, float, float,
                                                               float, float, float, float, float, float, float]):
        self.residue = residue
        self.acc_sym_values = acc_sym_values


class BFactorInfo:
    def __init__(self, residue: Residue, b_ave: float, b_c_alpha: float, bbb: float, b_side: float, b_outer: float):
        self.residue = residue
        self.b_ave = b_ave
        self.b_c_alpha = b_c_alpha
        self.bbb = bbb
        self.b_side = b_side
        self.b_outer = b_outer


class MoleculeValue:
    def __init__(self, begin_residue: Residue, end_residue: Residue, score: float):
        self.begin_residue = begin_residue
        self.end_residue = end_residue
        self.score = score


class ProfileEntry:
    def __init__(self, residue_type: str, profile_value: float):
        self.residue_type = residue_type
        self.profile_value = profile_value


class ResiduePlusProfile:
    def __init__(self, residue: Residue, profile: List[ProfileEntry], reliability: float):
        self.residue = residue
        self.profile = profile
        self.reliability = reliability


class ResiduePlusValue:
    def __init__(self, residue: Residue, value: float):
        self.residue = residue
        self.value = value


class ResiduePlusCount:
    def __init__(self, residue: Residue, count: int):
        self.residue = residue
        self.count = count


class ResiduePlustValuePlusAtoms:
    def __init__(self, residue: Residue, value: float, atoms: List[Atom]):
        self.residue = residue
        self.value = value
        self.atoms = atoms


class ResiduePlusValuePlusResidues:
    def __init__(self, residue: Residue, value: float, residues: List[Residue]):
        self.residue = residue
        self.value = value
        self.residues = residues


class SequenceInfo:
    def __init__(self, residue: Residue):
        self.residue = residue


class DSSPvalue:
    def __init__(self, residue: Residue, value: float):
        self.residue = residue
        self.value = value


class FileDBData:
    def __init__(self,
                 n_models: int,
                 resolution: float,
                 n_residues: int,
                 n_amino_acids: int,
                 n_h2o: int,
                 n_ligand_atoms: int,
                 pct_residues_with_bad_atoms: float,
                 pct_amino_acids_with_bad_backbone: float,
                 ramachandran_zscore: float,
                 chi1_chi2_zscore: float,
                 packing_quality: float,
                 rmsz_bonds: float,
                 rmsz_angles: float,
                 inside_outside_distribution: float,
                 bump_score: float,
                 use_this_file: bool):

        self.n_models = n_models
        self.resolution = resolution
        self.n_residues = n_residues
        self.n_amino_acids = n_amino_acids
        self.n_h2o = n_h2o
        self.n_ligand_atoms = n_ligand_atoms
        self.pct_residues_with_bad_atoms = pct_residues_with_bad_atoms
        self.pct_amino_acids_with_bad_backbone = pct_amino_acids_with_bad_backbone
        self.ramachandran_zscore = ramachandran_zscore
        self.chi1_chi2_zscore = chi1_chi2_zscore
        self.packing_quality = packing_quality
        self.rmsz_bonds = rmsz_bonds
        self.rmsz_angles = rmsz_angles
        self.inside_outside_distribution = inside_outside_distribution
        self.bump_score = bump_score
        self.use_this_file = use_this_file


# ==============================================================
# ==                                                          ==
# == Reusable PROCEDUREs                                      ==
# ==                                                          ==
# ==============================================================


def read_residue(line: str) -> Residue:
    fields = line.split(';')

    d = {
        "number": int(fields[0]),
        "pdb_number": int(fields[2]),
        "chain": fields[4],
        "type": fields[1],
    }

    if fields[3] != '':
        d["insertion_code"] = fields[3]

    if fields[5] != '_':
        d['model_number'] = int(fields[5])

    return d


def read_molecule_value(line: str) -> MoleculeValue:
    fields = line.split('|')
    return {"begin_residue": read_residue(fields[0]),
            "end_residue": read_residue(fields[1]),
            "value": float(fields[2])}


def read_profile_entry(line: str) -> ProfileEntry:
    fields = line.split(';')
    return {"amino_acid": fields[0], "value": float(fields[1])}


def read_atom(line: str) -> Atom:
    fields = line.split(';')

    d = {"residue": read_residue(line), "type": fields[6]}

    if fields[7] != "":
        d["alternative_atom"] = fields[7]

    return d


def read_coordinate(line: str) -> Coordinate:
    fields = line.split(',')
    return {"x": float(fields[0]), "y": float(fields[1]), "z": float(fields[2])}


def read_dot(line: str) -> Dot:
    return {"coordinate": {"x": float(fields[0]), "y": float(fields[1]), "z": float(fields[2])},
            "hue": int(fields[3])}


def read_line(line: str) -> Line:
    fields = line.split(';')

    return {"from": {"x": float(fields[0]), "y": float(fields[1]), "z": float(fields[2])},
            "to": {"x": float(fields[3]), "y": float(fields[4]), "z": float(fields[5])},
            "hue": int(fields[3])}


# ==============================================================
# ==                                                          ==
# == API Callables                                            ==
# ==                                                          ==
# ==============================================================

@bp.route('/test_empty/<pdbid>/', methods=['GET'])
def test_empty(pdbid: str):
    "Returns a quick nil-response"

    run_whatif(pdbid, "WEMPTY", False)

    return ""


@bp.route('/upload_pdb/', methods=['POST'])
def upload_pdb() -> str:
    """
    To use WHAT-IF functions on your own PDB files, you have to upload
    them first. The result of this function is an ID that you can use in
    further calls to wiwsd services.
    """

    pdbid = str(uuid4())

    pdb_path = os.path.join(flask.current_app.config['PDB_UPLOAD_DIRECTORY_PATH'], f"{pdbid}.pdb")

    flask.request.files['pdb'].save(pdb_path)

    return pdbid


@bp.route('/download_pdb/<pdbid>/', methods=['GET'])
def download_pdb(pdbid: str) -> str:
    """
    You can use this function to download PDB files. If you use the regular
    four letter code, you can retrieve standard PDB files, if you use the ID as
    returned by UploadPDB you can download a previously uploaded PDB file.
    """

    if len(pdbid) == 4:
        pdb_path = os.path.join(flask.current_app.config['PDB_DIRECTORY_PATH'], f"pdb{pdbid.lower()}.ent")
    else:
        pdb_path = os.path.join(flask.current_app.config['PDB_UPLOAD_DIRECTORY_PATH'], f"{pdbid}.pdb")

    with open(pdb_path, 'rt') as f:
        return f.read()


@bp.route('/get_surface_dots/<pdbid>/', methods=['GET'])
def get_surface_dots(pdbid: str) -> str:
    """
    Calculates the positions of dots that are homogeneously distributed over the surface
    of the molecule (waters are not incorporated in the calculation).
    Output is an array of structures that contain arrays filled with dots. A dot consist
    of three coordinates (x,y,z) and a number from 0-360 indicating the HUE value of the colour of
    the dot.
    Colours are on the Hue-wheel:
       0 = Blue
      60 = Purple
     120 = Red
     150 = Orange
     180 = Yellow
     240 = Green
     300 = Light Green
     360 = Blue (is zero again)
    """

    surface_dots = []
    dots = []
    last_residue_number = None
    for line in run_whatif("WSVDOT", pdbid, False):
        fields = line.split('|')
        atom = read_atom(fields[0])
        dot = read_dot(fields[1])
        if atom.residue.number != last_residue_number:

            if last_residue_number is not None:
                surface_dots.append({"atom": atom, "dots": dots})
                dots = []

            last_residue_number = atom.residue.number

        dots.append(dot)

    return flask.jsonify({"surface_dots": surface_dots})


@bp.route('/get_atom_accessibility_solvent/<pdbid>/', methods=['GET'])
def get_atom_accessibility_solvent(pdbid: str) -> str:
    """
    Returns for each atom in the input file its solvent accessibility
    in Å².
    Waters are neglected by this service.
    This service returns the 'accessible surface'.

    The atom_accessibility_molecular that runs the
    WHAT IF option WSVACM, will calculate the 'accessible molecular surface'.
    """

    infos = []
    for line in run_whatif("WSVACC", pdbid, False):
        fields = line.split('|')
        infos.append({"atom": read_atom(fields[0]), "accessibililty": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_atom_accessibility_molecular/<pdbid>/', methods=['GET'])
def get_atom_accessibility_molecular(pdbid: str) -> str:
    """
    Returns for each atom in the input file its accessibile molecular surface in Å².
    Waters are neglected by this service.
    This service returns the 'accessible molecular surface'. The atom_accessibility_molecular service, that runs the
    WHAT IF option WSVACC, will calculate the 'accessible surface'.
    The computation method is described in:
    """

    infos = []
    for line in run_whatif("WSVACM", pdbid, False):
        fields = line.split('|')
        infos.append({"atom": read_atom(fields[0]), "accessibililty": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_atom_accessibility_solvent_plus/<pdbid>/', methods=['GET'])
def get_atom_accessibility_solvent_plus(pdbid: str) -> str:
    """
    Returns for each atom in the input file its solvent accessibility in  Å².
    Waters are neglected by this service. A boolean is added to each atom indicating whether this atom
    is considered 'correct' by WHAT IF or not. Occupancy=0.0 is the most common reason for this
    boolean flag being returned as false.
    This service returns the 'accessible surface'. The atom_accessibility_solvent_plus service that runs the
    WHAT IF option WSVAMP, will calculate the 'accessible molecular surface'.
    """

    infos = []
    for line in run_whatif("WSVACP", pdbid, False):
        fields = line.split('|')
        infos.append({"atom": read_atom(fields[0]), "accessibililty": float(fields[1]), "status": fields[2] == "T"})

    return flask.jsonify({"infos": infos})


@bp.route('/get_atom_accessibility_molecular_plus/<pdbid>/', methods=['GET'])
def get_atom_accessibility_molecular_plus(pdbid: str) -> str:
    """
    Returns for each atom in the input file its accessibile molecular surface in Å²
    Waters are neglected by this service. A boolean is added to each atom indicating whether this atom
    is considered 'correct' by WHAT IF or not. Occupancy=0.0 is the most common reason for this
    boolean flag being returned as false.
    This service returns the 'accessible molecular surface'. The atom_accessibility_molecular service, that runs the
    WHAT IF option WSVACP, will calculate the 'accessible surface'.
    """

    infos = []
    for line in run_whatif("WSVAMP", pdbid, False):
        fields = line.split('|')
        infos.append({"atom": read_atom(fields[0]), "accessibililty": float(fields[1]), "status": fields[2] == "T"})

    return flask.jsonify({"infos": infos})


@bp.route('/get_residue_accessibility_molecular/<pdbid>/', methods=['GET'])
def get_residue_accessibility_molecular(pdbid: str) -> str:
    """
    Returns for each residue the backbone solvent accessibility, sidechain solvent accessibility,
    total solvent accessibility (sum of other two), plus a flag that is T or F. T when
    all atoms in this residue are OK, F when one or more are not OK.
    This service returns the 'accessible molecular surface'. The residue_accessibility_solvent service, that runs
    the WHAT IF option WSVAR2, will do the same but with 'accessible surface'.
    """

    infos = []
    for line in run_whatif("WSVAM2", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "backbone_accessibililty": float(fields[1]),
                      "side_chain_accessibililty": float(fields[2]),
                      "total_accessibililty": float(fields[3]),
                      "status": fields[4] == "T"})

    return flask.jsonify({"infos": infos})


@bp.route('/get_total_accessibility_molecular/<pdbid>/', methods=['GET'])
def get_total_accessibility_molecular(pdbid: str) -> str:
    """
    Returns for each entity the backbone, side chain, and total
    solvent accessibility, plus a flag that is T or F. T when
    all atoms in this entity are OK, F when one or more are not OK.
    For sugars and ligands the backbone is not defined, so their backbone accessibility is zero, and the
    side chain and total accessibility are the same.
    This service returns the &quot;accessible molecular surface&quot;. The tot_accessibility_solvent service, that runs
    the WHAT IF option TSVAR2, will do the same but with 'accessible surface'.
    The accessible surface is always larger than the accessible molecular surface.
    """

    infos = []
    for line in run_whatif("TSVAM2", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "accessibility": float(fields[1]),
                      "status": fields[2] == "T"})

    return flask.jsonify({"infos": infos})


@bp.route('/get_residue_accessibility_solvent/<pdbid>/', methods=['GET'])
def get_residue_accessibility_solvent(pdbid: str) -> str:
    """
    Returns for each residue the backbone molecular accessible surface, sidechain molecular accessible surface,
    total molecular accessible surface (sum of other two), plus a flag that is T or F. T when
    all atoms in this residue are OK, F when one or more are not OK.
    This service returns the 'accessbile surface'. The residue_accessibilty_solvent service, that runs
    the WHAT IF option WSVAM2, will do the same but with 'accessbile molecular surface'.
    """

    infos = []
    for line in run_whatif("WSVAR2", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "backbone_accessibililty": float(fields[1]),
                      "side_chain_accessibililty": float(fields[2]),
                      "total_accessibililty": float(fields[3]),
                      "status": fields[4] == "T"})

    return flask.jsonify({"infos": infos})


@bp.route('/get_total_accessibility_solvent/<pdbid>/', methods=['GET'])
def get_total_accessibility_solvent(pdbid: str) -> str:
    """
    Returns for each entity the backbone, side chain, and total
    molecular accessibility, plus a flag that is T or F. T when
    all atoms in this entity are OK, F when one or more are not OK.
    For sugars and ligands the backbone is not defined, so their backbone accessibility is zero, and the
    side chain and total accessibility are the same.
    This service returns the 'molecular surface'. The tot_accessibility_solvent service, that runs
    the WHAT IF option TSVAM2, will do the same but with 'solvent accessible surface'.
    """

    infos = []
    for line in run_whatif("TSVAR2", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "accessibililty": float(fields[1]),
                      "status": fields[2] == "T"})

    return flask.jsonify({"infos": infos})


@bp.route('/get_residue_accessibility_vacuum/<pdbid>/', methods=['GET'])
def get_residue_accessibility_vacuum(pdbid: str) -> str:
    """
    Returns for each residue the accessible surface, and the vacuum accessible surface. The latter
    is defined as the accessibility of the residue when taken out of the protein together with the
    backbone atoms of any residue it is covalently bound to.
    This service returns the 'accessible surface'. the residue_accessibility_vacuum_molecular service, that runs
    the WHAT IF option WSVAM3, will do the same but with 'accessible molecular surface'.
    The accessible surface is always larger than the accessible molecular surface.
    """

    infos = []
    for line in run_whatif("WSVAR3", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "normal_accessibility": float(fields[1]),
                      "vacuum_accessibility": float(fields[2])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_residue_accessibility_vacuum_molecular/<pdbid>/', methods=['GET'])
def get_residue_accessibility_vacuum_molecular(pdbid: str) -> str:
    """
    Returns for each residue the molecular surface, and the vacuum molecular surface. The latter
    is defined as the molecular surface of the residue when taken out of the protein together with the
    backbone atoms of any residue it is covalently bound to.
    This service returns the 'molecular surface'. The ResidueAccessibilityVacuum service, that runs
    the WHAT IF option WSVAR3, will do the same but with 'accessible surface'.
    """

    infos = []
    for line in run_whatif("WSVAM3", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "normal_accessibility": float(fields[1]),
                      "vacuum_accessibility": float(fields[2])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_residue_torsions_backbone/<pdbid>/', methods=['GET'])
def get_residue_torsions_backbone(pdbid: str) -> str:
    """
    Returns for each residue all its backbone torsion angles: phi, psi, omega, in degrees.
    They are returned as 999 if they don't exist, like phi that doesn't exist in an N-terminal residue.
    """

    infos = []
    for line in run_whatif("WSVCHB", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "phi_angle": float(fields[1]),
                      "psi_angle": float(fields[2]),
                      "omega_angle": float(fields[3])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_residue_torsions/<pdbid>/', methods=['GET'])
def get_residue_torsions(pdbid: str) -> str:
    """
    Returns for each residue all its torsion angles: phi, psi, omega, chi in degrees.
    They are returned as 999 if they don't exist, chi-angles are listed only if they
    exist and are returned as 999.9 when they are not meaningful (i.e. when one of the atoms defining
    defining that chi-angle is missing from the PDB file or has all coordinates zero.
    """

    infos = []
    for line in run_whatif("WSVCHI", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "phi_angle": float(fields[1]),
                      "psi_angle": float(fields[2]),
                      "omega_angle": float(fields[3]),
                      "chi_angles": [float(fields[4]), float(fields[5]), float(fields[6]), float(fields[7])]})

    return flask.jsonify({"infos": infos})


@bp.route('/get_cysteine_torsions/<pdbid>/', methods=['GET'])
def get_cysteine_torsions(pdbid: str) -> str:
    """
    Returns for each cysteine (bridge) all its torsion angles (in degrees).
    These are the torsion angles over the following five bonds, respectively:
    1. C-alpha - C-beta
    2. C-beta - S-gamma
    3. S-gamma - S-gamma
    4. S-gamma - C-beta
    5. C-beta - C-alpha
    """

    infos = []
    for line in run_whatif("WSVCYT", pdbid, False):
        fields = line.split('|')
        infos.append({"residue1": read_residue(fields[0]),
                      "residue2": read_residue(fields[1]),
                      "alpha1_beta1_angle": float(fields[2]),
                      "beta1_gamma1_angle": float(fields[3]),
                      "gamma1_gamma2_angle": float(fields[4]),
                      "gamma2_beta2_angle": float(fields[5]),
                      "beta2_alpha2_angle": float(fields[6])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_cysteine_bridges/<pdbid>/', methods=['GET'])
def get_cysteine_bridges(pdbid: str) -> str:
    """
    Show Cysteine bridges.
    """

    infos = []
    for line in run_whatif("WSVCYS", pdbid, False):
        fields = line.split('|')
        infos.append({"residue1": read_residue(fields[0]),
                      "residue2": read_residue(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_cysteine_bridges_no_pisa/<pdbid>/', methods=['GET'])
def get_cysteine_bridges_no_pisa(pdbid: str) -> str:
    """
    Show Cysteine bridges, but skip those that are between Pisa
    generated copies.
    This option was needed for one particular HOPE application...
    """

    infos = []
    for line in run_whatif("WSVCNP", pdbid, False):
        fields = line.split('|')
        infos.append({"residue1": read_residue(fields[0]),
                      "residue2": read_residue(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_cysteines_free/<pdbid>/', methods=['GET'])
def get_cysteines_free(pdbid: str) -> str:
    """
    Show Cysteine bridges.
    A cysteine is called free when it neither is involved in a cysteine bridge,
    nor functions as a ligand to a metal.
    """

    infos = []
    for line in run_whatif("WSVCY1", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_cysteines_metal/<pdbid>/', methods=['GET'])
def get_cysteines_metal(pdbid: str) -> str:
    """
    Show cysteines that are bound to a metal.
    Metals are either included when they are covalently bound to a cysteine (as
    can sometimes be observed for Pb, for example), or when the cysteine S-gamma is a ligand
    of the metal. The output consists of a row of cysteine-metal pairs.
    """

    infos = []
    for line in run_whatif("WSVCYM", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]), "metal": read_residue(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_bumps/<pdbid>/', methods=['GET'])
def get_bumps(pdbid: str) -> str:
    """
    Show Interatomic bumps
    The output are those pairs of atoms considered a bump, and their degree of interpenetration.
    A bump is listed if the Van der Waals' radii of two atoms interpenetrate more
    than 0.25 Å
    """

    infos = []
    for line in run_whatif("WSVBMP", pdbid, False):
        fields = line.split('|')
        infos.append({"residue1": read_residue(fields[0]),
                      "residue2": read_residue(fields[1]),
                      "bump_severity": float(fields[2])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_contacts_normal/<pdbid>/', methods=['GET'])
def get_contacts_normal(pdbid: str) -> str:
    """
    For each amino acid its atomic contacts are counted with everything but
    water. All distances are scored using a way too simple force-field.
    Two atoms are in contact when there is less than 0.5 Å
    between their Van der Waals surfaces.
    These contact scores are added up for all side chain atoms, and
    returned to the caller as a single value that relates (poorly) to the
    effects if the side chain studied is mutated.
    """

    infos = []
    for line in run_whatif("WSLCON", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_residues = int(fields[2])

        residues = []
        for i in range(3, n_residues + 3):
            residues.append(read_residue(fields[i]))

        infos.append({"residue": read_residue(fields[0]),
                      "value": value,
                      "residues": residues})

    return flask.jsonify({"infos": infos})


@bp.route('/get_contacts_relaxed/<pdbid>/', methods=['GET'])
def get_contacts_relaxed(pdbid: str) -> str:
    """
    For each amino acid its atomic contacts are counted with everything but
    water. All distances are scored using a way too simple force-field.
    Two atoms are in contact when there is less than 1.5 Å
    between their Van der Waals surfaces.
    These contact scores are added up for all  side chain atoms, and
    returned to the caller as a single value that relates (poorly) to the
    effects if the side chain studied is mutated.
    """

    infos = []
    for line in run_whatif("WSLCNR", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_residues = int(fields[2])

        residues = []
        for i in range(3, n_residues + 3):
            residues.append(read_residue(fields[i]))

        infos.append({"residue": read_residue(fields[0]),
                      "value": value,
                      "residues": residues})

    return flask.jsonify({"infos": infos})


@bp.route('/get_side_chain_contacts_normal/<pdbid>/', methods=['GET'])
def get_side_chain_contacts_normal(pdbid: str) -> str:
    """
    For each amino acid its side chain contacts (excluding C-beta)
    are counted with everything but water. All distances are scored
    using a way too simple force-field.
    Two atoms are in contact when there is less than 0.5 Å
    between their Van der Waals surfaces.
    These contact scores are added up for all side chain atoms and
    returned to the caller as a single value that relates (poorly) to the
    effects if the side chain studied is mutated.
    """

    infos = []
    for line in run_whatif("WSLCN2", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_residues = int(fields[2])

        residues = []
        for i in range(3, n_residues + 3):
            residues.append(read_residue(fields[i]))


        infos.append({"residue": read_residue(fields[0]),
                      "value": value,
                      "residues": residues})

    return flask.jsonify({"infos": infos})


@bp.route('/get_side_chain_contacts_relaxed/<pdbid>/', methods=['GET'])
def get_side_chain_contacts_relaxed(pdbid: str) -> str:
    """
    For each amino acid its side chain contacts (excluding C&beta;) are counted with everything but water.
    All distances are scored using a way too simple force-field.
    Two atoms are in contact when there is less than 1.5 Å between their Van der Waals surfaces.
    These contact scores are added up
    for all side chain atoms, and returned to the caller as a single value
    that relates (poorly) to the effects if the side chain studied is mutated.
    """

    infos = []
    for line in run_whatif("WSLCN3", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_residues = int(fields[2])

        residues = []
        for i in range(3, n_residues + 3):
            residues.append(read_residue(fields[i]))

        infos.append({"residue": read_residue(fields[0]),
                      "value": value,
                      "residues": residues})

    return flask.jsonify({"infos": infos})


@bp.route('/get_multimer_contacts/<pdbid>/', methods=['GET'])
def get_multimer_contacts(pdbid: str) -> str:
    """
    For each amino acid in the input file the contacts with
    other protein chains in the multimer are counted.
    This is a simple count, and has no linear relation to
    any form of contact energy
    """

    infos = []
    for line in run_whatif("WSVANP", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_nucleic_contacts/<pdbid>/', methods=['GET'])
def get_nucleic_contacts(pdbid: str) -> str:
    """
    For each amino acid in the input file the contacts with nucleic
    acids are counted.acids are counted.
    This is a simple count, and has no linear relation to any form of
    contact energycontact energy
    """

    infos = []
    for line in run_whatif("WSVHDC", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_protein_nucleic_contacts/<pdbid>/', methods=['GET'])
def get_protein_nucleic_contacts(pdbid: str) -> str:
    """
    Show protein - nucleic acid interactions
    The output are those pairs of atoms of which one is protein and the other nucleic acid, that
    have less than 1.0 Å between their Van der Waals' surfaces.
    """

    infos = []
    for line in run_whatif("WSVANU", pdbid, False):
        fields = line.split('|')

        infos.append({"atom1": read_atom(fields[0]),
                      "atom2": read_atom(fields[1]),
                      "contact_distance": float(fields[2])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_metal_contacts/<pdbid>/', methods=['GET'])
def get_metal_contacts(pdbid: str) -> str:
    """
    For each amino acid in the input file the ion contacts are counted.
    Depending on the distance they get a score ranging from 1.0 for close
    contacts to 0.25 for weak contacts. These contact scores are added up
    for all ions for all side chain atoms, and returned to the caller as a single value
    that relates (poorly) to the energetic contribution of this particular
    contact to the total folding energy.
    """

    infos = []
    for line in run_whatif("WSVHIC", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_metal_contacts_plus/<pdbid>/', methods=['GET'])
def get_metal_contacts_plus(pdbid: str) -> str:
    """
    For each amino acid in the input file the ion contacts are counted.
    Depending on the distance they get a score ranging from 1.0 for close
    contacts to 0.25 for weak contacts. These contact scores are added up
    for all ions for all side chain atoms, and returned to the caller as a single value
    that relates (poorly) to the energetic contribution of this particular
    contact to the total folding energy.
    The contacted metal ions are listed too.
    """

    infos = []
    for line in run_whatif("WSLHIC", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_atom = int(fields[2])

        atoms = []
        for i in range(3, n_atom + 3):
            atoms.append(read_atom(fields[i]))

        infos.append({"residue": residue,
                      "value": value,
                      "atoms": atoms})

    return flask.jsonify({"infos": infos})


@bp.route('/get_negative_ion_contacts/<pdbid>/', methods=['GET'])
def get_negative_ion_contacts(pdbid: str) -> str:
    """
    For each amino acid in the input file the ion contacts are counted.
    Depending on the distance they get a score ranging from 1.0 for close
    contacts to 0.25 for weak contacts. These contact scores are added up
    for all ions for all side chain atoms, and returned to the caller as a single value
    that relates (poorly) to the energetic contribution of this particular
    contact to the total folding energy.
    """

    infos = []
    for line in run_whatif("WSVHIN", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_negative_ion_contacts_plus/<pdbid>/', methods=['GET'])
def get_negative_ion_contacts_plus(pdbid: str) -> str:
    """
    For each amino acid in the input file the ion contacts are counted.
    Depending on the distance they get a score ranging from 1.0 for close
    contacts to 0.25 for weak contacts. These contact scores are added up
    for all ions for all side chain atoms, and returned to the caller as a single value
    that relates (poorly) to the energetic contribution of this particular
    contact to the total folding energy.
    The contacted metal ions are listed too.
    """

    infos = []
    for line in run_whatif("WSLHIN", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_residues = int(fields[2])

        residues = []
        for i in range(3, n_residues + 3):
            residues.append(read_residue(fields[i]))

        infos.append({"residue": residue,
                      "value": value,
                      "residues": residues})

    return flask.jsonify({"infos": infos})


@bp.route('/get_drug_contacts/<pdbid>/', methods=['GET'])
def get_drug_contacts(pdbid: str) -> str:
    """
    Show drug-like ligand interactions
    The output are those pairs of atoms of which one is a drug like ligand and the other macromolecule, that
    have less than 1.0 Å between their Van der Waals' surfaces.
    Compared to the get_ligand_contacts (WSVLIC) service, this one essentially skips contacts with ions.
    """

    infos = []
    for line in run_whatif("WSVLC2", pdbid, False):
        fields = line.split('|')

        infos.append({"residue1": read_residue(fields[0]),
                      "residue2": read_residue(fields[1]),
                      "contact_distance": float(fields[2])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_drug_contacts_short/<pdbid>/', methods=['GET'])
def get_drug_contacts_short(pdbid: str) -> str:
    """
    Show drug-like ligand interactions
    The output lists all residues that have a contact with (at least one) a drug-like ligand, and
    the shortest atom-to-atom amino acid - ligand contact distance, and all ligands this amino acid makes a contact with.
    A residue and a ligand are called in contact if they share a pair of atoms that
    have less than 1.0 Å between their Van der Waals' surfaces.
    Compared to the get_ligand_contacts (WSVLIC) service, this one essentially skips contacts with
    ions, sugars, crystallisation additives, and other small things.
    """

    infos = []
    for line in run_whatif("WSVLC3", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_residues = int(fields[2])

        residues = []
        for i in range(3, n_residues + 3):
            residues.append(read_residue(fields[i]))

        infos.append({"residue": residue,
                      "value": value,
                      "residues": residues})

    return flask.jsonify({"infos": infos})


@bp.route('/get_ligand_contacts/<pdbid>/', methods=['GET'])
def get_ligand_contacts(pdbid: str) -> str:
    """
    The output are those pairs of atoms of which one is ligand and the other macromolecule, that
    have less than 1.0 &Aring;ngstr&ouml;m between their Van der Waals' surfaces.
    This service considers drugs, metabolites, lipids, co-factors, and ions all as ligands. ATP
    ADP, NADPH, etc are considered as ligands unless they have the same name as nucleotides.
    """

    infos = []
    for line in run_whatif("WSVLIC", pdbid, False):
        fields = line.split('|')

        infos.append({"atom1": read_atom(fields[0]),
                      "atom2": read_atom(fields[1]),
                      "contact_distance": float(fields[2])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_salt_bridges/<pdbid>/', methods=['GET'])
def get_salt_bridges_residues(pdbid: str) -> str:
    """
    For each amino acid in the input file its side chain salt bridges are counted, valued, and added up.
    The total quasi salt bridge energy is returned to the caller as a single value
    that probably relates poorly to the energetic contribution of this particular
    salt bridge to the total folding energy. So, a salt bridge between the N-terminus and a Glu
    side chain is counted for the Glu, but not for the N-terminal residue, because we only check
    side chain salt bridges.
    Two residues are said to have a saltbridge if a negative atom in the one residue is within
    8.0 Å of a positive atom in the other residue.
    """

    infos = []
    for line in run_whatif("WSVHSB", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_salt_bridges_plus/<pdbid>/', methods=['GET'])
def get_salt_bridges_residues_plus(pdbid: str) -> str:
    """
    For each amino acid in the input file its side chain salt bridges are counted, valued, and added up.
    The total quasi salt bridge energy is returned to the caller as a single value
    that probably relates poorly to the energetic contribution of this particular
    salt bridge to the total folding energy. So, a salt bridge between the N-terminus and a Glu
    side chain is counted for the Glu, but not for the C-terminal residue, because we only check
    side chain salt bridges. All saltbridge partners are listed too.
    Two residues are said to have a saltbridge if a negative atom in the one residue is within
    8.0 Å of a positive atom in the other residue.
    """

    infos = []
    for line in run_whatif("WSLHSB", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        value = float(fields[1])
        n_residues = int(fields[2])

        residues = []
        for i in range(3, n_residues + 3):
            residues.append(read_residue(fields[i]))

        infos.append({"residue": residue,
                      "value": value,
                      "residues": residues})

    return flask.jsonify({"infos": infos})


@bp.route('/get_hydrogen_bonds_m/<pdbid>/', methods=['GET'])
def get_hydrogen_bonds_m(pdbid: str) -> str:
    """
    Lists potential hydrogenbonds between amino acids and anything but water.
    See get_hydrogen_bonds if you want to include hydrogenbonds with water too.
    The hydrogen bonding network is not optimised, and Asn, Gln, and His flips
    are not worked out prior to analyzing the hydrogen bonds.
    Output is a list of records, each record consisting of:
    The atoms making the hydrogenbond
    Four hydrogenbond geometric parameters:
     1) Donor Acceptor distance
     2) Proton Acceptor distance
     3) The deviation from 180 degrees over the acceptor
     4) The deviation from 180 degrees over the proton
    In case of bi-directional hydrogen bonds, both are listed, one-after-the-other
    """

    infos = []
    for line in run_whatif("WSVHB-", pdbid, False):
        fields = line.split('|')

        atom1 = read_atom(fields[0])
        atom2 = read_atom(fields[1])
        hydrogen_bond_parameters = [float(value) for value in fields[2].split(';') if len(value) > 0]

        infos.append({"atom1": atom1, "atom2": atom2,
                      "hydrogen_bond_parameters": hydrogen_bond_parameters})

    return flask.jsonify({"infos": infos})


@bp.route('/get_hydrogen_bonds_side/<pdbid>/', methods=['GET'])
def get_hydrogen_bonds_side(pdbid: str) -> str:
    """
    For each amino acid in the input file the H-bond enthalpy is calculated
    for each side chain atom. All these enthalpies are added up and the total
    enthalpy is returned.
    Residues that have no hydrogen bonding side chain atoms
    get the value 0.0. The hydrogenbonding network is optimized (see references) before the
    hydrogenbonds are analysed. Surface located water molecules are deleted before the network optimisation
    to make this a CPU-technically doable web service. Contact us if you need this facility with
    water molecules included, perhaps we can work something out.
    """

    infos = []
    for line in run_whatif("WSVHBF", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_hydrogen_bonds_no_water/<pdbid>/', methods=['GET'])
def get_hydrogen_bonds_no_water(pdbid: str) -> str:
    """
    For each amino acid in the input file the H-bond enthalpy is calculated
    for each side chain atom. All these enthalpies are added up and the total
    enthalpy is returned.
    Residues that have no hydrogen-bonding side chain atoms
    get the value 0.0. The hydrogenbonding network is optimized (see references) before the
    hydrogenbonds are analysed. All water molecules are deleted before the network optimisation
    to make this a CPU-technically doable web service. Contact us if you need this facility with
    water molecules included, perhaps we can work something out.
    """

    infos = []
    for line in run_whatif("WSVHBN", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_salt_bridges/<pdbid>/', methods=['GET'])
def get_salt_bridges(pdbid: str) -> str:
    """
    The output are those pairs of atoms considered a salt-bridge,
    and their inter-atomic distance. Histidines are not included in this calculation.
    """

    infos = []
    for line in run_whatif("WSVSBR", pdbid, False):
        fields = line.split('|')

        infos.append({"atom1": read_atom(fields[0]),
                      "atom2": read_atom(fields[1]),
                      "distance": float(field[3])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_salt_bridges_h/<pdbid>/', methods=['GET'])
def get_salt_bridges_h(pdbid: str) -> str:
    """
    Show saltbridges
    The output are those pairs of atoms considered a salt-bridge,
    and their inter-atomic distance. All histidines are considered
    positively charged in this calculation.
    """

    infos = []
    for line in run_whatif("WSVSB2", pdbid, False):
        fields = line.split('|')

        infos.append({"atom1": read_atom(fields[0]),
                      "atom2": read_atom(fields[1]),
                      "distance": float(field[3])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_hydrogen_bonds/<pdbid>/', methods=['GET'])
def get_hydrogen_bonds(pdbid: str) -> str:
    """
    Lists potential hydrogenbonds between amino acids and anything.
    The hydrogen bonding network is not optimised, and Asn, Gln, and His flips
    are not worked out prior to analyzing the hydrogen bonds.
    Output is a list of records, each record consisting of:
    The atoms making the hydrogenbond
    Four hydrogenbond geometric parameters:
     1) Donor Acceptor distance
     2) Proton Acceptor distance
     3) The deviation from 180 degrees over the acceptor
     4) The deviation from 180 degrees over the proton
    In case of bi-directional hydrogen bonds, both are listed one-after-the-other
    """

    infos = []
    for line in run_whatif("WSVHBO", pdbid, False):
        fields = line.split('|')

        infos.append({"atom1": read_atom(fields[0]),
                      "atom2": read_atom(fields[1]),
                      "hydrogenbond_parameters": [float(value) for value in fields[2].split(';') if len(value) > 0]})

    return flask.jsonify({"infos": infos})


# ==============================================================
# ==                                                          ==
# == Structure quality options                                ==
# ==                                                          ==
# ==============================================================


@bp.route('/get_packing_quality/<pdbid>/', methods=['GET'])
def get_packing_quality(pdbid: str) -> str:
    """
    Returns for each residue is packing normality. The normal score for proteins solved
    correctly at high resolution is about -0.4. A score of -5.0 or lower indicates a (very) poorly
    packed residue. If low scores are seen more than a few times, or a series of times
    in a row, one should start worrying about those residues.
    This work is described in:
    Quality control of protein models: Directional atomic contact analysis.
    G. Vriend, C. Sander. J.Appl.Cryst. (1993) 26, 47-60.
    """

    infos = []
    for line in run_whatif("WSVQUA", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_packing_quality_molecule/<pdbid>/', methods=['GET'])
def get_packing_quality_molecule(pdbid: str) -> str:
    """
    Returns packing normality for whole entry. The normal score for proteins solved
    correctly at high resolution is about -0.4. A score of -3.0 or lower indicates a (very) poorly
    packed molecule.
    This work is described in:
    Quality control of protein models: Directional atomic contact analysis.
    G. Vriend, C. Sander. J.Appl.Cryst. (1993) 26, 47-60.
    """

    infos = []
    for line in run_whatif("WSVQUA", pdbid, False):
        infos.append(read_molecule_value(line))

    return flask.jsonify({"infos": infos})


@bp.route('/get_inside_outside_distribution/<pdbid>/', methods=['GET'])
def get_inside_outside_distribution(pdbid: str) -> str:
    """
    Returns an inside - outside normality score that is callibrated against PDB files of
    water-soluble molecules. One score is returned for each molecule > 50 amino acids and
    one score for the whole file. Per molecule(s) scored the first residue, the last residue,
    and a score are listed
    """

    infos = []
    for line in run_whatif("WSVINO", pdbid, False):
        infos.append(read_molecule_value(line))

    return flask.jsonify({"infos": infos})


@bp.route('/get_improper_quality_sum/<pdbid>/', methods=['GET'])
def get_improper_quality_sum(pdbid: str) -> str:
    """
    Returns for each residue the sum of all improper dihedral Z-scores observable
    in that residue. One thus expects larger residues to generally score higher
    than small residues.
    """

    infos = []
    for line in run_whatif("WSVCHN", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_improper_quality_max/<pdbid>/', methods=['GET'])
def get_improper_quality_max(pdbid: str) -> str:
    """
    Returns for each residue the maximum of all its improper dihedral Z-scores.
    """

    infos = []
    for line in run_whatif("WSVCHX", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


# ==============================================================
# ==                                                          ==
# == Typical X-ray options                                    ==
# ==                                                          ==
# ==============================================================

@bp.route('/get_residue_b_factors/<pdbid>/', methods=['GET'])
def get_residue_b_factors(pdbid: str) -> str:
    """
    Returns per amino acid (in this order) the values:

    Average B-factor of all atoms;
    B-factor of the alpha carbon;
    Average B-factor of four backbone atoms;
    Average B-factor of the side chain atoms;
    Average B-factor of the outer (maximally 4) side chain atoms;

    Please be aware that for glycine the C-alpha is used as side chain, for alanine
    the C-beta is its whole side chain; and for many (small) amino acid types will the latter
    two values be identical.
    Only atoms are used that seem intact. If no intact atoms are found for any of the categories
    the value of 99.99 will be returned.
    """

    infos = []
    for line in run_whatif("WSVBFA", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "b_factor_all": float(fields[1]),
                      "b_factor_alphas": float(fields[2]),
                      "b_factor_backbone": float(fields[3]),
                      "b_factor_side_chains": float(fields[4]),
                      "b_factor_outer": float(fields[5])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_tau_angle/<pdbid>/', methods=['GET'])
def get_tau_angle(pdbid: str) -> str:
    """
    For each canonical amino acid in the input file the backbone angle tau
    is listed.
    """

    infos = []
    for line in run_whatif("WSVTAU", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "tau_angle": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_symmetry_contact/<pdbid>/', methods=['GET'])
def get_symmetry_contact(pdbid: str) -> str:
    """
    Returns for each residue the number of symmetry contacts made by that
    whole residue.
    Two atoms (of which one in 'another' asymmetric unit) are called in
    contact if the distance between their Van der Waals' surfaces is less
    than 0.25 Å
    """

    infos = []
    for line in run_whatif("WSVSMC", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "contact_count": int(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_symmetry_contact_total/<pdbid>/', methods=['GET'])
def get_symmetry_contact_total(pdbid: str) -> str:
    """
    Returns for each thing the number of symmetry contacts made by that
    whole thing. A thing can be amino acid, nucleic acid, water, sugar,
    ligand, ion, drug, etc.
    Two atoms (of which one in 'another' asymmetric unit) are called in
    contact if the distance between their Van der Waals' surfaces is less
    than 0.25 Å
    """

    infos = []
    for line in run_whatif("TSVSMC", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "contact_count": int(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_symmetry_contact_side/<pdbid>/', methods=['GET'])
def get_symmetry_contact_side(pdbid: str) -> str:
    """
    Returns for each amino acid the number of symmetry contacts made by
    it's side chain. For glycine the C-alpha is considered its side chain.
    Two atoms (of which one in 'another' asymmetric unit) are called in
    contact if the distance between their Van der Waals' surfaces is less
    than 0.25 Å.
    Contacts with symmetry related water molecules are not taken into account.
    """

    infos = []
    for line in run_whatif("WSVSMS", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "contact_count": int(fields[1])})

    return flask.jsonify({"infos": infos})


# ==============================================================
# ==                                                          ==
# == Protein engineering options                              ==
# ==                                                          ==
# ==============================================================


@bp.route('/get_mutation_bumps/<pdbid>/', methods=['GET'])
def get_mutation_bumps(pdbid: str) -> str:
    """
    Returns for each residue in the mutated protein the integrated bump
    value after optimally executing the mutatation.
    Two atoms (of which one in 'another' asymmetric unit) are called
    bumping if their Van der Waals' surfaces penetrate more than
    0.5 Å.
    This Web service accepts as input information about a mutation
    described by a residue number as integer, i.e. the number associated
    to this residue in the PDB file, and the new residue type as a string.
    Residue types can be given in 1-letter code (ACDEFGHIKLMNPQRSTVWY) or
    in 3-letter code (Ala, Cys, Asp, etc). Residue types are not case
    sensitive.
    The bump calculation is preceded by the execution of the requested mutation.
    """

    infos = []
    for line in run_whatif("WSVMBA", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_profile_hssp/<pdbid>/', methods=['GET'])
def get_profile_hssp(pdbid: str) -> str:
    """
    Warning. This service only returns results for PDB files known by HSSP, and
    those are the protein-containing subset of what the wwPDB distributes.
    This option reads the Profile from HSSP files, and returns those
    profile values. It also returns a thing called reliability. This is simply the
    number of sequences aligned at this position divided by the total number of sequences in the HSSP alignment.
    """

    infos = []
    for line in run_whatif("WSVPRF", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])

        profile_entries = []
        for i in range(1, 21):
            profile_entries.append(read_profile_entry(fields[i]))

        infos.append({"residue": residue,
                      "profile": profile_entries,
                      "reliability": float(fields[21])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_variability_hssp/<pdbid>/', methods=['GET'])
def get_variability_hssp(pdbid: str) -> str:
    """
    Warning. This service only returns results for PDB files known by HSSP, and
    those are the protein-containing subset of what the wwPDB distributes.
    This option reads the VAR column from HSSP files, and returns those
    values
    """

    infos = []
    for line in run_whatif("WSVVAR", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_likely_rotamers/<pdbid>/', methods=['GET'])
def get_likely_rotamers(pdbid: str) -> str:
    """
    Show Rotamer likelyhoods for all 20 amino acid types at each position in the protein.
    Output is a list of records, each record consisting of:
     * the residue
     * likelihoods for the 20 amino acid types (ACDEFGHIKLMNPQRSTVWY)
     * an estimator for the reliability of the 20 likelihoods for this residue
    This web-service is rather time consuming and can only be used for 100 amino acids
    at a time. You can determine which range of residues you want WHAT IF to analyze with the
    from and to parameters. Be aware that you cannot choose a range larger than 100 amino acids;
    and if you choose a larger range, WHAT IF will simply truncate your range. If you give non-existent
    residue numbers, WHAT IF will return no results.
    """

    infos = []
    for line in run_whatif("WSVROT", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        rotamer_table = [float(value) for value in fields[1].split(';') if len(value) > 0]
        quality_estimate = float(fields[2])

        infos.append({"residue": residue,
                      "rotamer_table": rotamer_table,
                      "quality_estimate": quality_estimate})

    return flask.jsonify({"infos": infos})


@bp.route('/get_access_and_symmetry/<pdbid>/', methods=['GET'])
def get_access_and_symmetry(pdbid: str) -> str:
    """
    Lists accessibilities and symmetry contact details
    Show accessibilities and symmetry contact details, decomposed in classes.
    Output is a list of records, each record consisting of:
     * the residue
     * 15 values:
       1) VACUUM ACCESSIBILITY
       2) ACCESSIBILITY IN OWN MOLECULE
       3) ACCESSIBILITY IN OWN MOLECULE AND ALL DRUGS
       4) ACCESSIBILITY IN CONTEXT OF WHOLE PDB FILE
       5) ACCESSIBILITY WITH SYMMETRY ON
       6 - 10) AS 1-5 BUT WITH PROBE RADIUS OF 1.0 Å
       11) NUMBER OF CONTACTS WITH RESIDUES; NO SYMMETRY
       12) NUMBER OF CONTACTS WITH NON-RESIDUE; NO SYMMETRY
       13) ADDITIONAL NUMBER OF CONTACTS WITH RESIDUES DUE TO SYMMETRY
       14) ADDITIONAL NUMBER OF CONTACTS WITH NON-RESIDUE DUE TO SYMMETRY
       15) RESERVED
    Value called RESERVED are for the time being returned as zero.
    """

    infos = []
    for line in run_whatif("WSVJOP", pdbid, False):
        fields = line.split('|')

        residue = read_residue(fields[0])
        acc_sym_values = [float(s) for s in fields[1].split(';')]

        infos.append({"residue": residue,
                      "values": acc_sym_values})

    return flask.jsonify({"infos": infos})


@bp.route('/get_proline_mutation_value/<pdbid>/', methods=['GET'])
def get_proline_mutation_value(pdbid: str) -> str:
    """
    This service determines for each position the molecule that is further than two positions away
    from a terminus a score that relates to the chance that a proline, when introduced at this
    position, would increase the stability of the whole protein.
    This option runs over all amino acids in the whole PDB file. If your PDB files holds too many
    amino acids (too many is somewhere between 500 and 1000) this option will time out and
    you should use proline_mutation_value_range several times instead.
    """

    infos = []
    for line in run_whatif("WSVPRO", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_happy_value/<pdbid>/', methods=['GET'])
def get_happy_value(pdbid: str) -> str:
    """
    This service determines for each position the molecule that is further than two positions away
    from a terminus a score that relates to the chance that a mutation at this position is
    could increase the stability of the whole protein.
    This option can time out when run over more than 500 amino acids. In those cases
    you should use this option several times instead giving the full range at once.
    """

    infos = []
    for line in run_whatif("WSVHAP", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


@bp.route('/get_proline_mutation_value_range/<pdbid>/<from_residue_number>/<to_residue_number>/', methods=['GET'])
def get_proline_mutation_value_range(pdbid: str, from_residue_number: int, to_residue_number: int) -> str:
    """
    This service determines for each position the given range (that is further than two positions away
    from any terminus) a score that relates to the chance that a proline, when introduced at this
    position, would increase the stability of the whole protein.
    """

    infos = []
    for line in run_whatif(f"WSVPRR {from_residue_number} {to_residue_number}", pdbid, False):
        fields = line.split('|')

        infos.append({"residue": read_residue(fields[0]),
                      "value": float(fields[1])})

    return flask.jsonify({"infos": infos})


# ==============================================================
# ==                                                          ==
# == Administrative options                                   ==
# ==                                                          ==
# ==============================================================

@bp.route('/get_pdb_sequence/<pdbid>/', methods=['GET'])
def get_pdb_sequence(pdbid: str) -> str:
    """
    The sequence in a PDB file is not uniquely defined. This server returns what
    WHAT IF thinks is the sequence in the PDB file you uploaded.
    """

    residues = []
    for line in run_whatif("WSVSEQ", pdbid, False):
        fields = line.split('|')

        residues.append(read_residue(fields[0]))

    return flask.jsonify({"residues": residues})


@bp.route('/get_het_group_names/<pdbid>/', methods=['GET'])
def get_het_group_names(pdbid: str) -> str:
    """
    The PDB has, for obscure reasons, decided to have atoms and hetatoms. Ligands, lipids, and
    ions tend to become hetatoms, and so do sometimes also groups that are attached to
    amino acids.
    The PDB tries to give every molecule that consists of hetatoms a unique three letter code
    and a unique name, but dont rely on this...
    This web-service returns for each molecule that consists of hetatoms its unique three
    letter code and its unique name. If a common name is given (a HETSYM name in
    PDB terms), this common name is returned rather than the formal name
    """

    infos = []
    for line in run_whatif("WSVHET", pdbid, False):
        fields = line.split('|')
        infos.append((fields[0], fields[1]))

    return flask.jsonify({"infos": infos})


@bp.route('/get_residues_dssp/<pdbid>/', methods=['GET'])
def get_residues_dssp(pdbid: str) -> str:
    """
    Returns for each residue the DSSP determined secondary structure in three-state (HSC).
    """

    infos = []
    for line in run_whatif("WSVHST", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "dssp_value": fields[1]})

    return flask.jsonify({"infos": infos})


@bp.route('/get_residues_dssp_plus_local/<pdbid>/', methods=['GET'])
def get_residues_dssp_plus_local(pdbid: str) -> str:
    """
    Returns for each residue the DSSP determined secondary structure in three-state (HSC).
    It than determines (coarsly) what the secondary structure roughly is according to
    the local phi-psi angles and adds a character (UHST). U for unknown or none of the others;
    H for helix; S for strand; T for the left handed helix area where phi and psi both are
    positive (or psi perhaps a bit negative but not much).
    """

    infos = []
    for line in run_whatif("WSVHSL", pdbid, False):
        fields = line.split('|')
        infos.append({"residue": read_residue(fields[0]),
                      "dssp_value": fields[1]})

    return flask.jsonify({"infos": infos})
