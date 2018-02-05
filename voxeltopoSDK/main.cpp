#include "synth_data_generator.h"
#include <reader.h> // tao volume reader

#include <iostream>
#include <string>
#include <map>
#include <experimental/filesystem>
#include <gflags/gflags.h>

using std::cout;
using std::endl;
using std::string;
namespace fs = std::experimental::filesystem;

/*all cmd line arguments*/
enum ValidMode {
	INTENSITY, SDF, NOISE
};
std::map<string, ValidMode> validmodes_map = {
	{"i", ValidMode::INTENSITY},
	{"s", ValidMode::SDF},
	{ "n", ValidMode::NOISE }
};
// define cmd flag 'md' for modes & its handler
DEFINE_string( md, "not-specified", "specify a mode for the utility to perform certain task. REQUIRED." );
static bool processMode( const char* _flagname, const std::string& _value )
{
	if ( validmodes_map.find( _value ) == validmodes_map.end() )
	{
		cout << "Invalid value for -md: " << _value << endl;
		return false;
	}
	return true;
}
DEFINE_validator( md, processMode );
// options for mode "i" and "s"
DEFINE_string( I, "not-specified", "input mesh for which certain function will be generated. REQUIRED." );
DEFINE_double( r, 0, "resolution of output volume. REQUIRED." );
// options for mode "n"
DEFINE_string( V, "not-specified", "input volume for which certain operation will be performed. REQUIRED." );
DEFINE_double( ns, 0.001, "noise scale (w.r.t. value range in given volume). OPTIONAL." );

void printUsage();
int main( int _argc, char ** _argv )
{
	google::SetUsageMessage( "\
		progexe -md=<a valid mode> <args...> \n\
		args depend on each mode. \n\
		" );
	google::ParseCommandLineFlags( &_argc, &_argv, true );

	if ( FLAGS_md == "i" || FLAGS_md == "s" )
	{
		/* figure out input file and construct volume's output filename */
		int res = FLAGS_r;
		fs::path mesh_filename = FLAGS_I;
		auto vol_filename = mesh_filename;
		vol_filename.replace_extension( ".mrc" );
		/* read in mesh */
		TriMesh *mesh = TriMesh::read( mesh_filename.string() );
		if ( !mesh )
		{
			cout << "Error: Cannot read mesh from file " << mesh_filename << endl;
			exit( -1 );
		}
		mesh->need_bbox();
		mesh->need_normals();
		synthdata::Generator::enforceNormalsOutward( mesh );

		/* now make volume */
		Volume vol( res, res, res );
		cout << "Volume space alloc.ed." << endl;
		if ( FLAGS_md == "i" )
		{
			cout << "computing space transform: volume -> mesh..." << endl;
			auto vol_box = synthdata::Generator::mkBBoxForVol( &vol );
			auto vol2mesh = synthdata::Generator::computeSpaceTransform( vol_box, mesh->bbox );
			cout << "making intensity volume" << endl;
			synthdata::Generator::mkIntensityVol( mesh, vol2mesh, &vol, synthdata::MVCInten() );
		}
		else if ( FLAGS_md == "s" )
		{
			cout << "making SDF volume" << endl;
			synthdata::Generator::mkSDFVol( mesh, &vol );
		}

		/* clean-up */
		delete mesh;

		/* output volume */
		cout << "writing volume to file: " << vol_filename << endl;
		int len = vol_filename.string().size() + 100;
		char* vol_filename_cpy = new char[ len ];
		memset( vol_filename_cpy, '\0', len );
		memcpy( vol_filename_cpy, vol_filename.string().c_str(), vol_filename.string().size() );
		vol.toMRCFile( vol_filename_cpy );
		delete vol_filename_cpy;
	}
	else if ( FLAGS_md == "n" )
	{
		fs::path vol_in_fn = FLAGS_V;
		fs::path vol_out_fn = vol_in_fn;
		auto ext = vol_in_fn.extension().string();
		vol_out_fn.replace_filename( vol_in_fn.filename().stem().string() + "_n" + ext );

		int len = vol_in_fn.string().size() + 100;
		char * tmp_fn = new char[ len ];
		memset( tmp_fn, '\0', len );
		memcpy( tmp_fn, vol_in_fn.string().c_str(), vol_in_fn.string().size() );
		MRCReader mrc_reader( tmp_fn );

		Volume* vol = mrc_reader.getVolume();
		double g_v = ( vol->getMax() - vol->getMin() ) * FLAGS_ns;
		printf( "adding noise (scale=%f) to volume... \n", g_v );
		synthdata::Generator::addNoise( vol, g_v );
		printf( "Done adding noise. \n" );
		len = vol_out_fn.string().size() + 100;
		memset( tmp_fn, '\0', len );
		memcpy( tmp_fn, vol_out_fn.string().c_str(), vol_out_fn.string().size() );
		printf( "Outputting volume to %s ...\n", tmp_fn );
		vol->toMRCFile( tmp_fn );
		printf( "Done outputting volume. \n" );
		delete tmp_fn;
		delete vol;
	}

	return 0;
}

void printUsage( int _argc, char ** _argv )
{
	cout << google::ProgramUsage() << endl;
}