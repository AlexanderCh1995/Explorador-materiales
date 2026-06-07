import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import py3Dmol

# ============================================================================
# CONFIGURACIÓN DE LA PÁGINA
# ============================================================================
st.set_page_config(
    page_title="Explorador de Materiales",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# BASE DE DATOS DE MATERIALES (extraída del Handbook de Shackelford)
# ============================================================================
MATERIALS_DB = {
    # METALES - Estructuras FCC
    "Aluminio (Al) - FCC": {
        "symbol": "Al", "name": "Aluminio", "atomic_number": 13,
        "atomic_mass": 26.982, "atomic_radius_nm": 0.143,
        "density_g_cm3": 2.70, "melting_point_C": 660.45,
        "crystal_system": "Cúbica Centrada en las Caras (FCC)",
        "lattice_params": {"a": 4.049, "b": 4.049, "c": 4.049, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Al", "Al", "Al", "Al"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        "color": "#C0C0C0", "category": "Metal",
        "description": "Metal ligero, excelente conductor, muy usado en aeronáutica.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/3/3d/Aluminium-4.jpg/320px-Aluminium-4.jpg"
    },
    "Cobre (Cu) - FCC": {
        "symbol": "Cu", "name": "Cobre", "atomic_number": 29,
        "atomic_mass": 63.546, "atomic_radius_nm": 0.128,
        "density_g_cm3": 8.93, "melting_point_C": 1084.87,
        "crystal_system": "Cúbica Centrada en las Caras (FCC)",
        "lattice_params": {"a": 3.615, "b": 3.615, "c": 3.615, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Cu", "Cu", "Cu", "Cu"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        "color": "#B87333", "category": "Metal",
        "description": "Excelente conductor eléctrico y térmico.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Copper_crystal.jpg/320px-Copper_crystal.jpg"
    },
    "Oro (Au) - FCC": {
        "symbol": "Au", "name": "Oro", "atomic_number": 79,
        "atomic_mass": 196.967, "atomic_radius_nm": 0.144,
        "density_g_cm3": 19.28, "melting_point_C": 1064.43,
        "crystal_system": "Cúbica Centrada en las Caras (FCC)",
        "lattice_params": {"a": 4.078, "b": 4.078, "c": 4.078, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Au", "Au", "Au", "Au"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        "color": "#FFD700", "category": "Metal",
        "description": "Metal precioso, altamente resistente a la corrosión.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/7/7e/Gold-crystals.jpg/320px-Gold-crystals.jpg"
    },
    "Plata (Ag) - FCC": {
        "symbol": "Ag", "name": "Plata", "atomic_number": 47,
        "atomic_mass": 107.868, "atomic_radius_nm": 0.144,
        "density_g_cm3": 10.50, "melting_point_C": 961.93,
        "crystal_system": "Cúbica Centrada en las Caras (FCC)",
        "lattice_params": {"a": 4.086, "b": 4.086, "c": 4.086, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Ag", "Ag", "Ag", "Ag"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        "color": "#E8E8E8", "category": "Metal",
        "description": "Mejor conductor eléctrico de todos los metales.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/9/9d/Silver_crystal.jpg/320px-Silver_crystal.jpg"
    },
    "Níquel (Ni) - FCC": {
        "symbol": "Ni", "name": "Níquel", "atomic_number": 28,
        "atomic_mass": 58.693, "atomic_radius_nm": 0.125,
        "density_g_cm3": 8.91, "melting_point_C": 1455,
        "crystal_system": "Cúbica Centrada en las Caras (FCC)",
        "lattice_params": {"a": 3.524, "b": 3.524, "c": 3.524, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Ni", "Ni", "Ni", "Ni"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        "color": "#B0B0B0", "category": "Metal",
        "description": "Metal resistente a la corrosión, usado en aleaciones.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/4/4d/Nickel_chunk.jpg/320px-Nickel_chunk.jpg"
    },
    "Plomo (Pb) - FCC": {
        "symbol": "Pb", "name": "Plomo", "atomic_number": 82,
        "atomic_mass": 207.2, "atomic_radius_nm": 0.175,
        "density_g_cm3": 11.34, "melting_point_C": 327.5,
        "crystal_system": "Cúbica Centrada en las Caras (FCC)",
        "lattice_params": {"a": 4.950, "b": 4.950, "c": 4.950, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Pb", "Pb", "Pb", "Pb"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        "color": "#5C5C5C", "category": "Metal",
        "description": "Metal denso y blando, usado en baterías y blindajes.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Lead_electrolytic_and_1cm3_cube.jpg/320px-Lead_electrolytic_and_1cm3_cube.jpg"
    },
    # METALES - Estructuras BCC
    "Hierro α (Fe) - BCC": {
        "symbol": "Fe", "name": "Hierro (α, BCC)", "atomic_number": 26,
        "atomic_mass": 55.845, "atomic_radius_nm": 0.124,
        "density_g_cm3": 7.87, "melting_point_C": 1538,
        "crystal_system": "Cúbica Centrada en el Cuerpo (BCC)",
        "lattice_params": {"a": 2.866, "b": 2.866, "c": 2.866, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Fe", "Fe"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "color": "#8A8A8A", "category": "Metal",
        "description": "Base del acero, material estructural fundamental.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Iron_electrolytic_and_1cm3_cube.jpg/320px-Iron_electrolytic_and_1cm3_cube.jpg"
    },
    "Cromo (Cr) - BCC": {
        "symbol": "Cr", "name": "Cromo", "atomic_number": 24,
        "atomic_mass": 51.996, "atomic_radius_nm": 0.125,
        "density_g_cm3": 7.19, "melting_point_C": 1863,
        "crystal_system": "Cúbica Centrada en el Cuerpo (BCC)",
        "lattice_params": {"a": 2.884, "b": 2.884, "c": 2.884, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Cr", "Cr"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "color": "#A8A8A8", "category": "Metal",
        "description": "Metal duro y brillante, usado en aceros inoxidables.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Chromium_crystals_and_1cm3_cube.jpg/320px-Chromium_crystals_and_1cm3_cube.jpg"
    },
    "Wolframio (W) - BCC": {
        "symbol": "W", "name": "Wolframio", "atomic_number": 74,
        "atomic_mass": 183.84, "atomic_radius_nm": 0.137,
        "density_g_cm3": 19.25, "melting_point_C": 3422,
        "crystal_system": "Cúbica Centrada en el Cuerpo (BCC)",
        "lattice_params": {"a": 3.165, "b": 3.165, "c": 3.165, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["W", "W"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "color": "#707070", "category": "Metal",
        "description": "Metal con el punto de fusión más alto de todos.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Tungsten_crystaline_structure.jpg/320px-Tungsten_crystaline_structure.jpg"
    },
    "Vanadio (V) - BCC": {
        "symbol": "V", "name": "Vanadio", "atomic_number": 23,
        "atomic_mass": 50.942, "atomic_radius_nm": 0.132,
        "density_g_cm3": 6.09, "melting_point_C": 1910,
        "crystal_system": "Cúbica Centrada en el Cuerpo (BCC)",
        "lattice_params": {"a": 3.024, "b": 3.024, "c": 3.024, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["V", "V"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "color": "#909090", "category": "Metal",
        "description": "Metal de transición, usado en aceros de alta resistencia.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Vanadium_crystals.jpg/320px-Vanadium_crystals.jpg"
    },
    # METALES - Estructuras HCP
    "Magnesio (Mg) - HCP": {
        "symbol": "Mg", "name": "Magnesio", "atomic_number": 12,
        "atomic_mass": 24.305, "atomic_radius_nm": 0.160,
        "density_g_cm3": 1.74, "melting_point_C": 650,
        "crystal_system": "Hexagonal Compacta (HCP)",
        "lattice_params": {"a": 3.209, "b": 3.209, "c": 5.211, "alpha": 90, "beta": 90, "gamma": 120},
        "species": ["Mg", "Mg"],
        "coords": [[0, 0, 0], [0.667, 0.333, 0.5]],
        "color": "#E8E8E8", "category": "Metal",
        "description": "Metal muy ligero, usado en aleaciones aeronáuticas.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Magnesium_crystals.jpg/320px-Magnesium_crystals.jpg"
    },
    "Zinc (Zn) - HCP": {
        "symbol": "Zn", "name": "Zinc", "atomic_number": 30,
        "atomic_mass": 65.38, "atomic_radius_nm": 0.133,
        "density_g_cm3": 7.13, "melting_point_C": 419.58,
        "crystal_system": "Hexagonal Compacta (HCP)",
        "lattice_params": {"a": 2.665, "b": 2.665, "c": 4.947, "alpha": 90, "beta": 90, "gamma": 120},
        "species": ["Zn", "Zn"],
        "coords": [[0, 0, 0], [0.667, 0.333, 0.5]],
        "color": "#B0B0B0", "category": "Metal",
        "description": "Usado en galvanización y aleaciones como el latón.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Zinc_crystal.jpg/320px-Zinc_crystal.jpg"
    },
    "Titanio (Ti) - HCP": {
        "symbol": "Ti", "name": "Titanio", "atomic_number": 22,
        "atomic_mass": 47.867, "atomic_radius_nm": 0.147,
        "density_g_cm3": 4.51, "melting_point_C": 1670,
        "crystal_system": "Hexagonal Compacta (HCP)",
        "lattice_params": {"a": 2.951, "b": 2.951, "c": 4.684, "alpha": 90, "beta": 90, "gamma": 120},
        "species": ["Ti", "Ti"],
        "coords": [[0, 0, 0], [0.667, 0.333, 0.5]],
        "color": "#C0C0C0", "category": "Metal",
        "description": "Metal fuerte y ligero, resistente a la corrosión.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Titanium_crystaline_structure.jpg/320px-Titanium_crystaline_structure.jpg"
    },
    # SEMICONDUCTORES - Estructura del Diamante
    "Silicio (Si) - Diamante": {
        "symbol": "Si", "name": "Silicio", "atomic_number": 14,
        "atomic_mass": 28.085, "atomic_radius_nm": 0.117,
        "density_g_cm3": 2.33, "melting_point_C": 1414,
        "crystal_system": "Cúbica del Diamante",
        "lattice_params": {"a": 5.431, "b": 5.431, "c": 5.431, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Si"] * 8,
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]],
        "color": "#3C3C3C", "category": "Semiconductor",
        "description": "Base de la electrónica moderna y los circuitos integrados.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Silicon_crystal.jpg/320px-Silicon_crystal.jpg"
    },
    "Germanio (Ge) - Diamante": {
        "symbol": "Ge", "name": "Germanio", "atomic_number": 32,
        "atomic_mass": 72.630, "atomic_radius_nm": 0.122,
        "density_g_cm3": 5.32, "melting_point_C": 938.3,
        "crystal_system": "Cúbica del Diamante",
        "lattice_params": {"a": 5.658, "b": 5.658, "c": 5.658, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Ge"] * 8,
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]],
        "color": "#505050", "category": "Semiconductor",
        "description": "Semiconductor usado en electrónica y fibra óptica.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Germanium_crystal.jpg/320px-Germanium_crystal.jpg"
    },
    # CERÁMICAS Y COMPUESTOS
    "Cloruro de Sodio (NaCl)": {
        "symbol": "NaCl", "name": "Cloruro de Sodio", "atomic_number": None,
        "atomic_mass": 58.44, "atomic_radius_nm": None,
        "density_g_cm3": 2.16, "melting_point_C": 801,
        "crystal_system": "Tipo NaCl (FCC)",
        "lattice_params": {"a": 5.640, "b": 5.640, "c": 5.640, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Na", "Na", "Na", "Na", "Cl", "Cl", "Cl", "Cl"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]],
        "color": "#FFFFFF", "category": "Cerámica",
        "description": "Sal común, estructura iónica típica.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/NaCl_crystal.jpg/320px-NaCl_crystal.jpg"
    },
    "Óxido de Magnesio (MgO)": {
        "symbol": "MgO", "name": "Óxido de Magnesio", "atomic_number": None,
        "atomic_mass": 40.304, "atomic_radius_nm": None,
        "density_g_cm3": 3.58, "melting_point_C": 2852,
        "crystal_system": "Tipo NaCl (FCC)",
        "lattice_params": {"a": 4.212, "b": 4.212, "c": 4.212, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Mg", "Mg", "Mg", "Mg", "O", "O", "O", "O"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]],
        "color": "#FFFFFF", "category": "Cerámica",
        "description": "Refractario de alta temperatura, estructura iónica.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Magnesium_oxide.jpg/320px-Magnesium_oxide.jpg"
    },
    "Fluorita (CaF₂)": {
        "symbol": "CaF₂", "name": "Fluorita", "atomic_number": None,
        "atomic_mass": 78.075, "atomic_radius_nm": None,
        "density_g_cm3": 3.18, "melting_point_C": 1418,
        "crystal_system": "Tipo Fluorita (FCC)",
        "lattice_params": {"a": 5.463, "b": 5.463, "c": 5.463, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Ca"]*4 + ["F"]*8,
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
                   [0.25, 0.25, 0.75], [0.75, 0.75, 0.75], [0.75, 0.25, 0.25], [0.25, 0.75, 0.25]],
        "color": "#90EE90", "category": "Cerámica",
        "description": "Mineral de fluoruro de calcio, usado en óptica.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Fluorite.jpg/320px-Fluorite.jpg"
    },
    "Cesio Cloruro (CsCl)": {
        "symbol": "CsCl", "name": "Cloruro de Cesio", "atomic_number": None,
        "atomic_mass": 168.36, "atomic_radius_nm": None,
        "density_g_cm3": 3.99, "melting_point_C": 645,
        "crystal_system": "Cúbica Simple",
        "lattice_params": {"a": 4.123, "b": 4.123, "c": 4.123, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Cs", "Cl"],
        "coords": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "color": "#FFFFFF", "category": "Cerámica",
        "description": "Estructura cúbica simple, coordinación 8:8.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Cesium_chloride.jpg/320px-Cesium_chloride.jpg"
    },
    "Zinc Blenda (ZnS)": {
        "symbol": "ZnS", "name": "Zinc Blenda (Esfalerita)",
        "atomic_number": None, "atomic_mass": 97.474,
        "atomic_radius_nm": None, "density_g_cm3": 4.09,
        "melting_point_C": 1850,
        "crystal_system": "Tipo Zinc Blenda (FCC)",
        "lattice_params": {"a": 5.406, "b": 5.406, "c": 5.406, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Zn"]*4 + ["S"]*4,
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]],
        "color": "#FFD700", "category": "Semiconductor",
        "description": "Estructura similar al diamante pero con dos elementos.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Sphalerite.jpg/320px-Sphalerite.jpg"
    },
    "Arseniuro de Galio (GaAs) - Zinc Blenda": {
        "symbol": "GaAs", "name": "Arseniuro de Galio", "atomic_number": None,
        "atomic_mass": 144.645, "atomic_radius_nm": None,
        "density_g_cm3": 5.32, "melting_point_C": 1238,
        "crystal_system": "Tipo Zinc Blenda (FCC)",
        "lattice_params": {"a": 5.653, "b": 5.653, "c": 5.653, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Ga"]*4 + ["As"]*4,
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]],
        "color": "#4B0082", "category": "Semiconductor",
        "description": "Semiconductor compuesto, usado en LEDs y láseres.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/GaAs_crystal.jpg/320px-GaAs_crystal.jpg"
    },
    "Carburo de Silicio (SiC)": {
        "symbol": "SiC", "name": "Carburo de Silicio", "atomic_number": None,
        "atomic_mass": 40.097, "atomic_radius_nm": None,
        "density_g_cm3": 3.21, "melting_point_C": 2730,
        "crystal_system": "Tipo Zinc Blenda (FCC)",
        "lattice_params": {"a": 4.360, "b": 4.360, "c": 4.360, "alpha": 90, "beta": 90, "gamma": 90},
        "species": ["Si"]*4 + ["C"]*4,
        "coords": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                   [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]],
        "color": "#2F4F4F", "category": "Cerámica",
        "description": "Material muy duro, usado en abrasivos y electrónica de potencia.",
        "image_url": "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0c/Silicon_carbide.jpg/320px-Silicon_carbide.jpg"
    },
}

# ============================================================================
# FUNCIONES AUXILIARES
# ============================================================================
@st.cache_data
def get_material_data(material_name):
    return MATERIALS_DB.get(material_name)

def create_3d_viewer(structure, width=500, height=400):
    cif_str = structure.to(fmt="cif")
    view = py3Dmol.view(width=width, height=height)
    view.addModel(cif_str, "cif")
    view.setStyle({"stick": {"radius": 0.08}, "sphere": {"scale": 0.3}})
    view.setBackgroundColor("white")
    view.zoomTo()
    return view._make_html()

# ============================================================================
# INTERFAZ DE USUARIO
# ============================================================================
st.title("🔬 Explorador de Materiales")
st.markdown("""
**Sistema de Información y Visualización de Materiales**
Basado en el *CRC Materials Science and Engineering Handbook* de Shackelford.
""")
st.divider()

st.sidebar.header("🔍 Selección de Material")
categories = sorted(set([m["category"] for m in MATERIALS_DB.values()]))
selected_category = st.sidebar.selectbox("Categoría:", ["Todos"] + categories)

if selected_category == "Todos":
    material_names = list(MATERIALS_DB.keys())
else:
    material_names = [name for name, data in MATERIALS_DB.items() if data["category"] == selected_category]

selected_material = st.sidebar.selectbox("Material:", material_names)

st.sidebar.divider()
st.sidebar.header("⚙️ Opciones")
show_xrd = st.sidebar.checkbox("Mostrar Difractograma XRD", value=True)
show_3d = st.sidebar.checkbox("Mostrar Estructura 3D", value=True)
xrd_min_angle = st.sidebar.slider("Ángulo 2θ mínimo (°)", 5, 30, 10)
xrd_max_angle = st.sidebar.slider("Ángulo 2θ máximo (°)", 60, 120, 90)

# ============================================================================
# CONTENIDO PRINCIPAL
# ============================================================================
if selected_material:
    data = get_material_data(selected_material)
    
    # CREACIÓN DE LA ESTRUCTURA CRISTALINA (Disponible para 3D, XRD y Densidad)
    lattice = Lattice.from_parameters(**data["lattice_params"])
    structure = Structure(lattice, data["species"], data["coords"])

    # --- SECCIÓN 1: INFORMACIÓN GENERAL ---
    st.header(f"📊 {data['name']} ({data['symbol']})")
    col1, col2, col3 = st.columns([1, 2, 1])

    with col1:
        st.subheader("🖼️ Imagen Representativa")
        try:
            st.image(data["image_url"], caption=data["name"], use_container_width=True)
        except:
            st.info("Imagen no disponible")
        st.markdown(f"**Categoría:** {data['category']}")

    with col2:
        st.subheader("📝 Descripción")
        st.markdown(data["description"])
        st.subheader("📋 Propiedades Fundamentales")
        props_data = {
            "Propiedad": ["Número Atómico (Z)", "Masa Atómica (amu)",
                         "Radio Atómico (nm)", "Densidad (g/cm³)",
                         "Punto de Fusión (°C)", "Sistema Cristalino"],
            "Valor": [
                str(data["atomic_number"]) if data["atomic_number"] else "N/A",
                f"{data['atomic_mass']:.3f}",
                f"{data['atomic_radius_nm']:.3f}" if data["atomic_radius_nm"] else "N/A",
                f"{data['density_g_cm3']:.2f}",
                f"{data['melting_point_C']:.1f}",
                data["crystal_system"]
            ]
        }
        st.table(pd.DataFrame(props_data))

    with col3:
        st.subheader("📐 Parámetros de Red")
        params = data["lattice_params"]
        params_data = {
            "Parámetro": ["a (Å)", "b (Å)", "c (Å)", "α (°)", "β (°)", "γ (°)"],
            "Valor": [f"{params['a']:.3f}", f"{params['b']:.3f}", f"{params['c']:.3f}",
                     f"{params['alpha']:.1f}", f"{params['beta']:.1f}", f"{params['gamma']:.1f}"]
        }
        st.table(pd.DataFrame(params_data))
        st.subheader("🔬 Información Cristalográfica")
        st.markdown(f"- **Átomos/iones por celda:** {len(data['species'])}")
        st.markdown(f"- **Especies:** {', '.join(set(data['species']))}")

    st.divider()

    # --- SECCIÓN 2: VISUALIZACIÓN 3D ---
    if show_3d:
        st.header("🧊 Estructura Cristalina 3D")
        col3d_1, col3d_2 = st.columns([2, 1])
        with col3d_1:
            html_3d = create_3d_viewer(structure)
            st.components.v1.html(html_3d, height=450)
            st.caption("🖱️ Usa el mouse para rotar, zoom y desplazar la estructura")
        with col3d_2:
            st.subheader("Información de la Celda")
            st.markdown(f"""
            - **Volumen de celda:** {structure.volume:.2f} Å³
            - **Número de átomos:** {len(structure)}
            - **Fórmula:** {structure.composition.reduced_formula}
            - **Densidad calculada:** {structure.density:.2f} g/cm³
            """)
            st.subheader("Posiciones Atómicas")
            positions = [{"Átomo": site.species_string,
                          "x": f"{site.frac_coords[0]:.3f}",
                          "y": f"{site.frac_coords[1]:.3f}",
                          "z": f"{site.frac_coords[2]:.3f}"} for site in structure]
            st.dataframe(pd.DataFrame(positions), use_container_width=True)

    st.divider()

    # --- SECCIÓN 3: DIFRACTOGRAMA XRD ---
    if show_xrd:
        st.header("📈 Difractograma de Rayos X (XRD)")
        st.markdown("""
        **Fundamento físico:** Ley de Bragg: $n\lambda = 2d\sin(\theta)$
        """)
        col_xrd_1, col_xrd_2 = st.columns([2, 1])
        with col_xrd_1:
            xrd_calc = XRDCalculator(wavelength=1.5406)
            pattern = xrd_calc.get_pattern(structure, two_theta_range=(xrd_min_angle, xrd_max_angle))
            fig = go.Figure()
            for i in range(len(pattern.x)):
                hkl = pattern.hkls[i]
                hkl_str = f"({hkl[0]} {hkl[1]} {hkl[2]})" if len(hkl) >= 3 else str(hkl)
                fig.add_trace(go.Scatter(
                    x=[pattern.x[i], pattern.x[i]], y=[0, pattern.y[i]],
                    mode='lines', line=dict(color='blue', width=2), name=hkl_str,
                    hoverinfo='text', text=f"2θ: {pattern.x[i]:.2f}°<br>Intensidad: {pattern.y[i]:.1f}<br>Plano: {hkl_str}"
                ))
            fig.update_layout(
                title=f"Patrón de Difracción - {data['name']}",
                xaxis_title="Ángulo 2θ (grados)", yaxis_title="Intensidad Relativa",
                xaxis=dict(range=[xrd_min_angle, xrd_max_angle]),
                yaxis=dict(range=[0, max(pattern.y) * 1.15]) if len(pattern.y) > 0 else None,
                template="plotly_white", showlegend=False, height=500
            )
            st.plotly_chart(fig, use_container_width=True)
        with col_xrd_2:
            st.subheader("Picos Principales")
            if len(pattern.x) > 0:
                peaks_data = {
                    "2θ (°)": [f"{x:.2f}" for x in pattern.x],
                    "d (Å)": [f"{1.5406/(2*np.sin(np.radians(x/2))):.3f}" for x in pattern.x],
                    "Intensidad": [f"{y:.1f}" for y in pattern.y],
                    "Plano (hkl)": [f"({h[0]} {h[1]} {h[2]})" if len(h)>=3 else str(h) for h in pattern.hkls]
                }
                st.dataframe(pd.DataFrame(peaks_data), use_container_width=True, hide_index=True)
            else:
                st.warning("No hay picos en este rango.")

    st.divider()

    # --- SECCIÓN 4: ANÁLISIS FÍSICO-MATEMÁTICO ---
    st.header("🧮 Análisis Físico-Matemático")
    col_analysis_1, col_analysis_2 = st.columns(2)

    with col_analysis_1:
        st.subheader("Factor de Empaquetamiento Atómico (APF)")
        if data["atomic_radius_nm"]:
            r = data["atomic_radius_nm"] * 10  # nm a Å
            n_atoms = len(data["species"])
            V_cell = structure.volume
            V_atoms = n_atoms * (4/3) * np.pi * r**3
            APF = V_atoms / V_cell
            st.markdown(f"""
            $$APF = \\frac{{n \\cdot \\frac{{4}}{{3}}\\pi r^3}}{{V_{{celda}}}}$$
            - **Radio atómico (r):** {r:.3f} Å
            - **Átomos por celda (n):** {n_atoms}
            - **Volumen de celda:** {V_cell:.2f} Å³
            - **APF calculado:** {APF:.3f}
            """)
            crystal_sys = data["crystal_system"]
            if "FCC" in crystal_sys: st.info(f"APF teórico FCC = 0.740 | Error: {abs(APF-0.740)*100:.1f}%")
            elif "BCC" in crystal_sys: st.info(f"APF teórico BCC = 0.680 | Error: {abs(APF-0.680)*100:.1f}%")
            elif "HCP" in crystal_sys: st.info(f"APF teórico HCP = 0.740 | Error: {abs(APF-0.740)*100:.1f}%")
            elif "Diamante" in crystal_sys: st.info(f"APF teórico Diamante = 0.340 | Error: {abs(APF-0.340)*100:.1f}%")
        else:
            st.warning("Radio atómico no disponible (compuestos iónicos).")

    with col_analysis_2:
        st.subheader("Densidad Teórica")
        # Mapa de masas atómicas
        atomic_masses = {
            "Al": 26.982, "Cu": 63.546, "Au": 196.967, "Ag": 107.868,
            "Ni": 58.693, "Pb": 207.2, "Fe": 55.845, "Cr": 51.996,
            "W": 183.84, "V": 50.942, "Mg": 24.305, "Zn": 65.38,
            "Ti": 47.867, "Si": 28.085, "Ge": 72.630, "Na": 22.990,
            "Cl": 35.453, "O": 15.999, "Ca": 40.078, "F": 18.998,
            "Cs": 132.905, "S": 32.065, "Ga": 69.723, "As": 74.922, "C": 12.011
        }
        
        try:
            # 1. Sumamos las masas atómicas de todos los átomos/iones dentro de la celda unitaria (n * A)
            total_mass = sum([atomic_masses[sp] for sp in data["species"]])
            
            # 2. Convertimos el volumen de Å³ a cm³ (1 Å = 10^-8 cm -> 1 Å³ = 10^-24 cm³)
            V_cell_cm3 = structure.volume * 1e-24 
            
            # 3. Número de Avogadro
            N_A = 6.022e23 
            
            # 4. Fórmula física correcta: ρ = (n * A) / (V_c * N_A)
            density_calc = total_mass / (V_cell_cm3 * N_A)

            st.markdown(f"""
            $$\\rho = \\frac{{n \\cdot A}}{{V_{{celda}} \\cdot N_A}}$$
            - **Masa total celda (n·A):** {total_mass:.2f} g/mol
            - **Volumen de celda (V_c):** {V_cell_cm3:.2e} cm³
            - **Densidad calculada:** {density_calc:.2f} g/cm³
            - **Densidad tabulada:** {data['density_g_cm3']:.2f} g/cm³
            - **Error relativo:** {abs(density_calc - data['density_g_cm3']) / data['density_g_cm3'] * 100:.2f}%
            """)
        except KeyError:
            st.warning("Faltan datos de masa atómica para algún elemento en el compuesto.")

    st.divider()

    # --- SECCIÓN 5: REFERENCIAS ---
    st.header("📚 Referencias")
    st.markdown("""
    - Shackelford, J.F. & Alexander, W. (2001). *CRC Materials Science and Engineering Handbook*. CRC Press.
    - Callister, W.D. & Rethwisch, D.G. (2014). *Materials Science and Engineering: An Introduction*. Wiley.
    - Datos cristalográficos: Materials Project Database (materialsproject.org)
    """)

st.divider()
st.markdown("""
<div style='text-align: center; color: gray;'>
    <p>🔬 Explorador de Materiales v2.1 | Desarrollado con Python, Streamlit, pymatgen y py3Dmol</p>
    <p>Basado en el CRC Materials Science and Engineering Handbook de Shackelford</p>
</div>
""", unsafe_allow_html=True)
