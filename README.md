##g_RotTransDNA
***

###About

This program calculates rotational/translational displacementof the DNA
around/along Z axis with respect to its starting structure (tpr).

***

###Requirements
To compile and install, GROMACS libraries <code> libgmx, libmd, libgmxana </code> are required.
***

###Download
<pre><code>git clone https://github.com/rjdkmr/g_RotTransDNA
</code></pre>
***

###Installation
<pre><code>cd g_RotTransDNA
mkdir build
cd build
cmake ..  -DGMX_PATH=/opt/gromacs
make
make install
</code></pre>

Directory <code>/opt/gromacs</code> should contains <code>include</code> and <code> lib </code> directories. If these directories are in seprate locations, use followings:
<pre><code>cmake ..  -DGMX_LIB=/path/to/lib -GMX_INCLUDE=/path/to/include -DCMAKE_INSTALL_PREFIX=/opt/gmxhole
</code></pre>

If fftw library <code> libfftw3f.so or libfftw3f.a </code> are not present in standard locations:
<pre><code>-DFFTW_LIB=/path/to/fftw3/lib</code></pre>
***

###Usage
<pre><code>g_RotTransDNA -h
</code></pre>
***
