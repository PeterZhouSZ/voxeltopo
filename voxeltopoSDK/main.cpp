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
	fs::path mesh_filename = fs::path( _argv[ 2 ] );
	auto vol_filename = mesh_filename; 
	vol_filename.replace_extension( ".mrc" );
	/* read in mesh */
	TriMesh *mesh = TriMesh::read( mesh_filename.string() );
	if ( !mesh )
	{
		cout << "Error: Cannot read mesh from file " << mesh_filename << endl;
		exit( -1 );
	}

	/* now make volume */
	Volume vol( res, res, res );
	cout << "Volume space alloc.ed." << endl;
	if ( mode.compare( "-i" ) == 0 )
	{
		cout << "making intensity volume" << endl;
		synthdata::Generator::mkIntensityVol( mesh, &vol );
	}
	else if ( mode.compare( "-s" ) == 0 )
	{
		cout << "making SDF volume" << endl;
		synthdata::Generator::mkSDFVol( mesh, &vol );
	}

	/* output volume */
	vol.toMRCFile( (char*)( vol_filename.c_str() ) );
	return 0;
}

void veryifyArgs( int _argc, char ** _argv )
{
	if ( _argc < 3 )
	{
		printUsage( _argc, _argv );
		exit( -1 );
	}
}

void printUsage( int _argc, char ** _argv )
{
	cout << "usage: " << _argv[ 0 ] << "<mode> <input-mesh-name>. " << endl 
		<< "mode can be - i( intensity ) or -s( sdf )." << endl;
}