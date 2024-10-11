from ifermi.interpolate import FourierInterpolator
from ifermi.kpoints import kpoints_from_bandstructure
from pymatgen.io.ase import AseAtomsAdaptor
from BoltzTraP2 import sphere
import numpy as np
from pymatgen.electronic_structure.core import Spin,OrbitalType
from pymatgen.symmetry.groups import SpaceGroup
from ifermi.surface import FermiSurface
from skimage import transform
import h5py    
from pymatgen.electronic_structure.plotter import BSDOSPlotter
import matplotlib.pyplot as plt
from plotly.offline import plot as show_plotly
from ifermi.plot import FermiSlicePlotter, FermiSurfacePlotter    
def write_superbands_atomate(h5pyfn,sb_id,formula,Tc):
    
    filename='superband/'+sb_id
    uniform_bs_entry = atomate_db.collection.find_one({'task_label': 'nscf uniform','formula_pretty': formula})
    bs_uniform = atomate_db.get_band_structure(uniform_bs_entry['task_id'])
 
    _kpoints = kpoints_from_bandstructure(bs_uniform)
    _atoms = AseAtomsAdaptor.get_atoms(bs_uniform.structure)
    equivalences = sphere.get_equivalences(atoms=_atoms,
                nkpt=_kpoints.shape[0] * 5, magmom=None)
    mesh = 2 * np.max(np.abs(np.vstack(equivalences)), axis=0) + 1
    # interpolate the energies onto a dense k-point mesh
    dense_bs, velocities = FourierInterpolator(bs_uniform).interpolate_bands(return_velocities=True)
    dense_kpoints = kpoints_from_bandstructure(dense_bs)
    _spins=dense_bs.bands.keys()
    if(len(_spins)==1):
        sc_bands=dense_bs.bands[Spin(1)]
    else:
        sc_bands=np.vstack((dense_bs.bands[Spin(1)],dense_bs.bands[Spin(-1)]))
    #取最中间12条能带
    lenth=18
    dense_bs_abs=np.zeros(len(sc_bands))
    if len(sc_bands)>lenth:
        for i in range(len(sc_bands)):
            dense_bs_abs[i]=sum(sc_bands[i]-dense_bs.efermi)
        argsorts=np.argsort(np.abs(dense_bs_abs))
        dense_bs_tmp=sc_bands[argsorts[:lenth]]-dense_bs.efermi
        dense_bs_abs=np.zeros(len(dense_bs_tmp))
        for i in range(len(dense_bs_tmp)):
            dense_bs_abs[i]=sum(dense_bs_tmp[i])
        argsorts=np.argsort(dense_bs_abs)

        #统一形状
        plot_z=dense_bs_tmp.reshape((-1,mesh[2],mesh[1],mesh[0]))
        bands_sc=np.zeros((lenth,32,32,32))
        for i in range(lenth):
            bands_sc[i]=transform.resize(plot_z[argsorts[i]],(32,32,32))
    else:
        bands_sc=np.zeros((lenth,32,32,32))
        startlen=(lenth-len(sc_bands))//2
        for i in range(len(sc_bands)):
            bands_sc[i+startlen]=transform.resize(sc_bands[i].reshape((mesh[2],mesh[1],mesh[0])),(32,32,32))

    SpaceGroups=SpaceGroup(bs_uniform.structure.get_space_group_info()[0])
    
    bs_uniform.structure.to(filename+'.cif')

    bs_dos = atomate_db.get_dos(uniform_bs_entry['task_id'])
    dos0s=bs_dos.energies-bs_dos.efermi
    for i in range(len(bs_dos.get_spd_dos().keys())):
        ors=bs_dos.get_spd_dos()[OrbitalType(i)]
        if len(ors.densities.keys())==1:
            dos0s=np.vstack((dos0s,ors.densities[Spin(1)]))
        else:
            dos0s=np.vstack((dos0s,(ors.densities[Spin(1)]+ors.densities[Spin(-1)])/2))
    try:        
        fs = FermiSurface.from_band_structure(
        dense_bs, mu=0.0, wigner_seitz=True, calculate_dimensionality=True,
        property_data=velocities, property_kpoints=dense_kpoints
        )

        fermi_lines_lens=[len(verts) for (verts, faces) in fs.all_vertices_faces()]
        verts=[vert for (vert, face) in fs.all_vertices_faces()]
        fermi_line=np.vstack(verts)
        
        #Fermi Furface 
        fig = plt.figure(figsize=(18, 7))
        ax1 = fig.add_subplot(121)
        FermiSlicePlotter(fs.get_fermi_slice((0, 0, 1),0.00)).get_plot(ax=ax1,cmin=100000,cmax=1000000)
        ax1.margins(y=0, x=0)
        ax2 = fig.add_subplot(122)
        FermiSlicePlotter(fs.get_fermi_slice((1, 0, 0),0.00)).get_plot(ax=ax2,cmin=100000,cmax=1000000)
        ax2.margins(y=0, x=0)
        fig.savefig(filename+"fs.pdf")

        plot = FermiSurfacePlotter(fs).get_plot(plot_type="plotly")
        show_plotly(plot, include_mathjax="cdn", 
                    filename=filename+'fs.html', auto_open=False,
                    include_plotlyjs='directory')
    except:
        print('FermiSurface Failure')
    
    try:
        line_bs_entry = atomate_db.collection.find_one({'task_label': 'nscf line', 'formula_pretty': formula})
        bs_sc = atomate_db.get_band_structure(line_bs_entry['task_id'])
        plt_1=BSDOSPlotter(bs_projection='elements', dos_projection='orbitals',vb_energy_range=5, fixed_cb_energy=5)
        plt_1.get_plot(bs=bs_sc,dos=bs_dos)
        plt.savefig(filename+"bs.pdf")
        #bs_sc.structure.to_file(filename+'.cif')
    except:
        print('Bandstructure Failure')    

    with h5py.File(h5pyfn, 'r+') as fast5_data:
        fast5_Read=fast5_data.create_group(sb_id)
        fast5_Read.attrs['formula']=_atoms.symbols.formula._formula
        fast5_Read.attrs['Tc']=Tc
        fast5_Read.attrs['SG_num']=SpaceGroups.int_number
        fast5_Read.attrs['SG_sys']=SpaceGroups.crystal_system
        fast5_Read.attrs['volume']=bs_uniform.structure.volume
        fast5_Read.create_dataset('atoms',data=_atoms.get_chemical_symbols())
        fast5_Read.create_dataset('cell',data=_atoms.cell)
        fast5_Read.create_dataset('position',data=_atoms.positions)
        fast5_Read.create_dataset('sc_bands',data=bands_sc)
        fast5_Read.create_dataset('sc_DOSs',data=dos0s.T)
        
        try:
            fast5_Read.create_dataset('fermi_line',data=fermi_line)
            fast5_Read.create_dataset('fermi_lens',data=fermi_lines_lens)
            fast5_Read.create_dataset('recip_latt',data=fs.reciprocal_space.reciprocal_lattice)
        except:
            print('Bandstructure Failure')    
    
    