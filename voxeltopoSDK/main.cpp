#include "synth_data_generator.h"
#include <iostream>
#include <string>
#include <experimental/filesystem>

using std::cout;
using std::endl;
using std::string;
namespace fs = std::experimental::filesystem;

void veryifyArgs( int _argc, char ** _argv );
void printUsage( int _argc, char ** _argv );
int main( int _argc, char ** _argv )
{
	veryifyArgs( _argc, _argv );
	/* figure out input file and construct volume's output filename */
	std::string mode = _argv[ 1 ];
	int res = std::stoi( _argv[ 2 ] );
	fs::path mesh_filename = fs::path( _argv[ 3 ] );
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
	
	/* now make volume */
	Volume vol( res, res, res );
	cout << "Volume space alloc.ed." << endl;
	if ( mode.compare( "-i" ) == 0 )
	{
		cout << "computing space transform: volume -> mesh..." << endl;
		auto vol_box = synthdata::Generator::mkBBoxForVol( &vol );
		auto vol2mesh = synthdata::Generator::computeSpaceTransform( vol_box, mesh->bbox );
		cout << "making intensity volume" << endl;
		synthdata::Generator::mkIntensityVol( mesh, vol2mesh, &vol, synthdata::SimpleInten() );
	}
	else if ( mode.compare( "-s" ) == 0 )
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

	return 0;
}

void veryifyArgs( int _argc, char ** _argv )
{
	if ( _argc < 4 )
	{
		printUsage( _argc, _argv );
		exit( -1 );
	}
}

void printUsage( int _argc, char ** _argv )
{
	cout << "usage: " << _argv[ 0 ] << "<mode> <vol-res> <input-mesh-name>. " << endl 
		<< "mode can be - i( intensity ) or -s( sdf )." << endl;
}