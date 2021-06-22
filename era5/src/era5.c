/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "driver.h"
#include "era5.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"
#include "netcdf.h"


#define alloc(p, s, t) {p = (t)malloc(sizeof(*p)*s);}
#define nc_catch(e) { \
    int e_ = e; \
    if (e_ != NC_NOERR) {\
        fprintf(stderr, "[%s: %d] %s\n", __FILE__, __LINE__, nc_strerror(e_)); \
        exit(EXIT_FAILURE); \
    }}
#ifdef SINGLE_PRECISION
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_float(a, b, c, d, e))
#else
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_double(a, b, c, d, e))
#endif


enum dimid
{
    TIME = 0,
    LEVEL,
    LAT,
    LON,
    LAYER,
    LW_WAVENUMBER,
    SW_WAVENUMBER,
    NUM_DIMS
};


struct Output
{
    int *dimid;
    int ncid;
    int *varid;
};

static char time_units[512];
static double *times;

static size_t global_num_lon;
static size_t nlon;
static double *lons;
static char lon_units[512];
static size_t nlat;
static double *lats;
static char lat_units[512];
static size_t nlevel;
static int x;
static int X;

/*Reorder from (t,z,y,x) to (t,y,x,z).*/
static void tzyx_to_tyxz(fp_t *dest, fp_t *src, int nx, int ny, int nz, int nt)
{
    int i;
    for (i=0; i<nt; ++i)
    {
        int j;
        for (j=0; j<nz; ++j)
        {
            int k;
            for (k=0; k<ny; ++k)
            {
                int m;
                for (m=0; m<nx; ++m)
                {
                    int offset1 = i*nz*ny*nx + j*ny*nx + k*nx + m;
                    int offset2 = i*ny*nx*nz + k*nx*nz + m*nz + j;
                    dest[offset2] = src[offset1];
                }
            }
        }
    }
    return;
}


/*Reorder from (t,y,x,z) to (t,z,y,x).*/
static void tyxz_to_tzyx(fp_t *dest, fp_t *src, int nx, int ny, int nz, int nt)
{
    int i;
    for (i=0; i<nt; ++i)
    {
        int j;
        for (j=0; j<ny; ++j)
        {
            int k;
            for (k=0; k<nx; ++k)
            {
                int m;
                for (m=0; m<nz; ++m)
                {
                    int offset1 = i*ny*nx*nz + j*nx*nz + k*nz + m;
                    int offset2 = i*nz*ny*nx + m*ny*nx + j*nx + k;
                    dest[offset2] = src[offset1];
                }
            }
        }
    }
    return;
}


/*Reserve memory and read in atmospheric data.*/
Atmosphere_t create_atmosphere(Parser_t *parser)
{
    /*Add/parse command line arguments.*/
    snprintf(parser->description, desclen,
             "Calculates the radiative fluxes for ERA5 reanalysis data.");
    add_argument(parser, "era5_file", NULL, "Input data file.", NULL);
    add_argument(parser, "ghg_file", NULL, "Greenhouse gase file.", NULL);
    int one = 1;
    add_argument(parser, "-CFC-11", NULL, "Path to CFC-11 cross sections.", &one);
    add_argument(parser, "-CFC-12", NULL, "Path to CFC-12 cross sections.", &one);
    add_argument(parser, "-CFC-113", NULL, "Path to CFC-113 cross sections.", &one);
    add_argument(parser, "-CH4", NULL, "Include methane.", NULL);
    add_argument(parser, "-clean", NULL, "Run without aerosols.", NULL);
    add_argument(parser, "-clear", NULL, "Run without clouds.", NULL);
    add_argument(parser, "-CO2", NULL, "Include carbon dioxide.", NULL);
    add_argument(parser, "-HCFC-22", NULL, "Path to HCFC-22 cross sections.", &one);
    add_argument(parser, "-H2O", NULL, "Include water vapor.", NULL);
    add_argument(parser, "-h2o-ctm", NULL, "Directory containing H2O continuum files", &one);
    add_argument(parser, "-N2-N2", NULL, "Path to N2-N2 CIA cross sections.", &one);
    add_argument(parser, "-N2O", NULL, "Include nitrous oxide.", NULL);
    add_argument(parser, "-O2-N2", NULL, "Path to O2-N2 CIA cross sections.", &one);
    add_argument(parser, "-O2-O2", NULL, "Path to O2-O2 CIA cross sections.", &one);
    add_argument(parser, "-O3", NULL, "Include O3.", NULL);
    add_argument(parser, "-o3-ctm", NULL, "Directory containing O3 continuum files", &one);
    add_argument(parser, "-t", "--time-lower-bound", "Starting time index.", &one);
    add_argument(parser, "-T", "--Time-upper-bound", "Ending time index.", &one);
    add_argument(parser, "-x", "--lon-lower-bound", "Starting longitude index.", &one);
    add_argument(parser, "-X", "--lon-upper-bound", "Ending longitude index.", &one);
    add_argument(parser, "-y", "--lat-lower-bound", "Starting latitude index.", &one);
    add_argument(parser, "-year", NULL, "Year for gas abundances.", &one);
    add_argument(parser, "-Y", "--lat-upper-bound", "Ending latitude index.", &one);
    add_argument(parser, "-z", "--level-lower-bound", "Starting level index.", &one);
    add_argument(parser, "-Z", "--level-upper-bound", "Ending level index.", &one);
    parse_args(*parser);

    /*Open the input file.*/
    char buffer[valuelen];
    get_argument(*parser, "era5_file", buffer);
    int ncid;
    nc_catch(nc_open(buffer, NC_NOWRITE, &ncid));

    /*Determine the number of times.*/
    Atmosphere_t atm;
    int t = get_argument(*parser, "-t", buffer) ? atoi(buffer) : 0;
    int T;
    if (get_argument(*parser, "-T", buffer))
    {
        T = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "time", &dimid));
        size_t num_times;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_times));
        T = (int)num_times - 1;
    }
    atm.num_times = T - t + 1;
    alloc(times, atm.num_times, double *);
    int varid;
    nc_catch(nc_inq_varid(ncid, "time", &varid));
    size_t start[4] = {t, 0, 0, 0};
    size_t count[4] = {atm.num_times, 1, 1, 1};
    get_var(ncid, varid, start, count, times);
    nc_catch(nc_get_att_text(ncid, varid, "units", time_units));

    /*Determine the number of columns.*/
    x = get_argument(*parser, "-x", buffer) ? atoi(buffer) : 0;
    if (get_argument(*parser, "-X", buffer))
    {
        X = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "lon", &dimid));
        size_t num_lon;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_lon));
        X = (int)num_lon - 1;
    }
    nlon = X - x + 1;
    alloc(lons, nlon, double *);
    nc_catch(nc_inq_varid(ncid, "lon", &varid));
    start[0] = x; start[1] = 0; start[2] = 0; start[3] = 0;
    count[0] = nlon; count[1] = 1; count[2] = 1; count[3] = 1;
    get_var(ncid, varid, start, count, lons);
    nc_catch(nc_get_att_text(ncid, varid, "units", lon_units));

    int y = get_argument(*parser, "-y", buffer) ? atoi(buffer) : 0;
    int Y;
    if (get_argument(*parser, "-Y", buffer))
    {
        Y = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "lat", &dimid));
        size_t num_lat;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_lat));
        Y = (int)num_lat - 1;
    }
    nlat = Y - y + 1;
    alloc(lats, nlat, double *);
    nc_catch(nc_inq_varid(ncid, "lat", &varid));
    start[0] = y; start[1] = 0; start[2] = 0; start[3] = 0;
    count[0] = nlat; count[1] = 1; count[2] = 1; count[3] = 1;
    get_var(ncid, varid, start, count, lats);
    nc_catch(nc_get_att_text(ncid, varid, "units", lat_units));
    atm.num_columns = nlon*nlat;

    /*Determine the number of levels.*/
    int z = get_argument(*parser, "-z", buffer) ? atoi(buffer) : 0;
    int Z;
    if (get_argument(*parser, "-Z", buffer))
    {
        Z = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "sigma_level", &dimid));
        size_t num_levels;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_levels));
        Z = (int)num_levels - 1;
    }
    nlevel = Z - z + 1;
    atm.num_levels = nlevel;
    atm.num_layers = atm.num_levels - 1;

    /*Pressure.*/
    alloc(atm.level_pressure, atm.num_times*atm.num_columns*atm.num_levels, fp_t *);
    fp_t *pressure;
    alloc(pressure, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    nc_catch(nc_inq_varid(ncid, "p", &varid));
    start[0] = t; start[1] = z; start[2] = y; start[3] = x;
    count[0] = atm.num_times; count[1] = atm.num_levels; count[2] = nlat; count[3] = nlon;
    get_var(ncid, varid, start, count, pressure);
    tzyx_to_tyxz(atm.level_pressure, pressure, nlon, nlat, atm.num_levels, atm.num_times);
    free(pressure);
    alloc(atm.layer_pressure, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
    int i;
    for (i=0; i<atm.num_times; ++i)
    {
        int j;
        for (j=0; j<atm.num_columns; ++j)
        {
            int k;
            for (k=0; k<atm.num_layers; ++k)
            {
                int offset1 = i*atm.num_columns*atm.num_layers + j*atm.num_layers + k;
                int offset2 = i*atm.num_columns*atm.num_levels + j*atm.num_levels + k;
                atm.layer_pressure[offset1] = 0.5*(atm.level_pressure[offset2] +
                                                   atm.level_pressure[offset2+1]);
            }
        }
    }

    /*Temperature.*/
    alloc(atm.level_temperature, atm.num_times*atm.num_columns*atm.num_levels, fp_t *);
    fp_t *temperature;
    alloc(temperature, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    nc_catch(nc_inq_varid(ncid, "t", &varid));
    start[0] = t; start[1] = z; start[2] = y; start[3] = x;
    count[0] = atm.num_times; count[1] = atm.num_levels; count[2] = nlat; count[3] = nlon;
    get_var(ncid, varid, start, count, temperature);
    tzyx_to_tyxz(atm.level_temperature, temperature, nlon, nlat, atm.num_levels, atm.num_times);
    free(temperature);
    alloc(atm.layer_temperature, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
    for (i=0; i<atm.num_times; ++i)
    {
        int j;
        for (j=0; j<atm.num_columns; ++j)
        {
            int k;
            for (k=0; k<atm.num_layers; ++k)
            {
                int layer_offset = i*atm.num_columns*atm.num_layers + j*atm.num_layers;
                int level_offset = i*atm.num_columns*atm.num_levels + j*atm.num_levels;
                fp_t *tlay = &(atm.layer_temperature[layer_offset]);
                fp_t const *tlev = &(atm.level_temperature[level_offset]);
                fp_t const *play = &(atm.layer_pressure[layer_offset]);
                fp_t const *plev = &(atm.level_pressure[level_offset]);
                tlay[k] = tlev[k] + (tlev[k+1] - tlev[k])*(play[k] - plev[k])/(plev[k+1] - plev[k]);
            }
        }
    }

    /*Molecular abundances.*/
    fp_t const to_ppmv = 1.e6;
    fp_t const dry_air_mass = 28.97;
    struct MoleculeMeta
    {
        int id;
        char *flag;
        char *name;
        fp_t mass;
    };
    int const num_molecules = 5;
    struct MoleculeMeta molecules[2] = {{H2O, "-H2O", "q", 18.01528},
                                        {O3, "-O3", "o3", 48.}};
    alloc(atm.molecules, num_molecules, int *);
    atm.num_molecules = 0;
    alloc(atm.ppmv, num_molecules, fp_t **);
    fp_t *abundance;
    alloc(abundance, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    for (i=0; i<2; ++i)
    {
        if (get_argument(*parser, molecules[i].flag, NULL))
        {
            atm.molecules[atm.num_molecules] = molecules[i].id;
            alloc(atm.ppmv[atm.num_molecules], atm.num_times*atm.num_columns*atm.num_levels, fp_t *);
            fp_t *ppmv = atm.ppmv[atm.num_molecules];
            nc_catch(nc_inq_varid(ncid, molecules[i].name, &varid));
            start[0] = t; start[1] = z; start[2] = y; start[3] = x;
            count[0] = atm.num_times; count[1] = atm.num_levels; count[2] = nlat; count[3] = nlon;
            get_var(ncid, varid, start, count, abundance);
            int j;
            for (j=0; j<atm.num_times*atm.num_levels*nlat*nlon; ++j)
            {
                abundance[j] *= to_ppmv*(dry_air_mass/molecules[i].mass);
            }
            tzyx_to_tyxz(ppmv, abundance, nlon, nlat, atm.num_levels, atm.num_times);
            atm.num_molecules++;
        }
    }
    free(abundance);

    /*Molecular continua.*/
    if (!get_argument(*parser, "-h2o-ctm", atm.h2o_ctm))
    {
        snprintf(atm.h2o_ctm, valuelen, "%s", "none");
    }
    if (!get_argument(*parser, "-o3-ctm", atm.o3_ctm))
    {
        snprintf(atm.o3_ctm, valuelen, "%s", "none");
    }
    atm.num_cfcs = 0;
    atm.num_cias = 0;
    atm.num_cia_species = 0;

    /*Surface temperature.*/
    alloc(atm.surface_temperature, atm.num_times*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "skt", &varid));
    start[0] = t; start[1] = y; start[2] = x; start[3] = 0;
    count[0] = atm.num_times; count[1] = nlat; count[2] = nlon; count[3] = 1;
    get_var(ncid, varid, start, count, atm.surface_temperature);

    /*Calculate the cosine of the solar zenith angle from the solar irradiance.*/
    int dimid;
    nc_catch(nc_inq_dimid(ncid, "lon", &dimid));
    nc_catch(nc_inq_dimlen(ncid, dimid, &global_num_lon));
    nc_catch(nc_inq_dimid(ncid, "lat", &dimid));
    size_t global_num_lat;
    nc_catch(nc_inq_dimlen(ncid, dimid, &global_num_lat));
    fp_t *lat;
    alloc(lat, global_num_lat, fp_t *);
    nc_catch(nc_inq_varid(ncid, "lat", &varid));
    start[0] = 0; start[1] = 0; start[2] = 0; start[3] = 0;
    count[0] = global_num_lat; count[1] = 1; count[2] = 1; count[3] = 1;
    get_var(ncid, varid, start, count, lat);
    fp_t *weights;
    alloc(weights, global_num_lat, fp_t *);
    fp_t total_weight = 0.;
    for (i=0; i<global_num_lat; ++i)
    {
        weights[i] = cos(2.*M_PI*lat[i]/360.);
        total_weight += weights[i];
    }
    free(lat);
    fp_t *irradiance;
    alloc(irradiance, atm.num_times*global_num_lat*global_num_lon, fp_t *);
    nc_catch(nc_inq_varid(ncid, "tisr", &varid));
    start[0] = t; start[1] = 0; start[2] = 0; start[3] = 0;
    count[0] = atm.num_times; count[1] = global_num_lat; count[2] = global_num_lon; count[3] = 1;
    get_var(ncid, varid, start, count, irradiance);
    fp_t *mean_irradiance;
    alloc(mean_irradiance, atm.num_times, fp_t *);
    fp_t const seconds_per_day = 86400.;
    for (i=0; i<atm.num_times; ++i)
    {
        mean_irradiance[i] = 0.;
        int j;
        for (j=0; j<global_num_lat; ++j)
        {
            fp_t running_sum = 0.;
            int k;
            for (k=0; k<global_num_lon; ++k)
            {
                int offset = i*global_num_lat*global_num_lon + j*global_num_lon + k;
                irradiance[offset] /= seconds_per_day;
                running_sum += irradiance[offset];
            }
            mean_irradiance[i] += (running_sum/global_num_lon)*(weights[j]/total_weight);
        }
        mean_irradiance[i] *= 4.;
    }
    free(weights);
    alloc(atm.solar_zenith_angle, atm.num_times*atm.num_columns, fp_t *);
    for (i=0; i<atm.num_times; ++i)
    {
        int j;
        for (j=0; j<nlat; ++j)
        {
            int k;
            for (k=0; k<nlon; ++k)
            {
                int offset1 = i*nlat*nlon + j*nlon + k;
                int offset2 = i*global_num_lat*global_num_lon + (j+y)*global_num_lon + k+x;
/*
                atm.solar_zenith_angle[offset1] = irradiance[offset2]/mean_irradiance[i];
*/
                atm.solar_zenith_angle[offset1] = -1.;
            }
        }
    }
    free(irradiance);
    free(mean_irradiance);

    /*Solar irradiance.*/
    alloc(atm.total_solar_irradiance, atm.num_times*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "tisr", &varid));
    start[0] = t; start[1] = y; start[2] = x; start[3] = 0;
    count[0] = atm.num_times; count[1] = nlat; count[2] = nlon; count[3] = 1;
    get_var(ncid, varid, start, count, atm.total_solar_irradiance);
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        atm.total_solar_irradiance[i] /= seconds_per_day*atm.solar_zenith_angle[i];
    }

    /*Surface albedo.*/
    atm.albedo_grid_size = 2;
    alloc(atm.albedo_grid, atm.albedo_grid_size, fp_t *);
    fp_t const ir_uv_boundary = 10000.;
    fp_t const ir_uv_offset = 1.e-5;
    atm.albedo_grid[0] = ir_uv_boundary - ir_uv_offset;
    atm.albedo_grid[1] = ir_uv_boundary + ir_uv_offset;
    alloc(atm.surface_albedo, atm.num_times*atm.num_columns*atm.albedo_grid_size, fp_t *);
    fp_t *albedo;
    alloc(albedo, atm.num_times*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "fal", &varid));
    start[0] = t; start[1] = y; start[2] = x; start[3] = 0;
    count[0] = atm.num_times; count[1] = nlat; count[2] = nlon; count[3] = 1;
    get_var(ncid, varid, start, count, albedo);
    for (i=0; i<atm.num_times*atm.num_columns; ++i)
    {
        atm.surface_albedo[i*atm.albedo_grid_size] = albedo[i];
        atm.surface_albedo[i*atm.albedo_grid_size + 1] = albedo[i];
    }
    free(albedo);


    /*Clouds.*/
    atm.clear = get_argument(*parser, "-clear", buffer) ? 1 : 0;
    if (!atm.clear)
    {
        /*Constants.*/
        double const gas_constant = 8.314462; /*[J mol-1 K-1].*/
        double const kg_per_g = 1./1000.; /*[kg g-1].*/
        double const molar_mass = 28.9647; /*[g mol-1].*/
        double const pa_per_mb = 100.; /*[Pa mb-1].*/

        /*Calclutate air density.*/
        fp_t *air_density;
        alloc(air_density, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
        for (i=0; i<atm.num_times*atm.num_columns*atm.num_layers; ++i)
        {
            air_density[i] = (atm.layer_pressure[i]*pa_per_mb*molar_mass*kg_per_g)/
                             (atm.layer_temperature[i]*gas_constant);
        }

        /*Read in the cloud properties.*/
        alloc(atm.cloud_fraction, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
        fp_t *clouds;
        alloc(clouds, atm.num_times*atm.num_layers*atm.num_columns, fp_t *);
        nc_catch(nc_inq_varid(ncid, "cc", &varid));
        start[0] = t; start[1] = z; start[2] = y; start[3] = x;
        count[0] = atm.num_times; count[1] = atm.num_layers; count[2] = nlat; count[3] = nlon;
        get_var(ncid, varid, start, count, clouds);
        tzyx_to_tyxz(atm.cloud_fraction, clouds, nlon, nlat, atm.num_layers, atm.num_times);
        for (i=0; i<atm.num_times*atm.num_columns*atm.num_layers; ++i)
        {
            if (atm.cloud_fraction[i] <= 0.)
            {
                atm.cloud_fraction[i] = 0.;
            }
        }
        alloc(atm.ice_water_content, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
        nc_catch(nc_inq_varid(ncid, "ciwc", &varid));
        start[0] = t; start[1] = z; start[2] = y; start[3] = x;
        count[0] = atm.num_times; count[1] = atm.num_layers; count[2] = nlat; count[3] = nlon;
        get_var(ncid, varid, start, count, clouds);
        tzyx_to_tyxz(atm.ice_water_content, clouds, nlon, nlat, atm.num_layers, atm.num_times);
        fp_t const g_per_kg = 1000.;
        for (i=0; i<atm.num_times*atm.num_columns*atm.num_layers; ++i)
        {
            if (atm.ice_water_content[i] <= 0.)
            {
                atm.ice_water_content[i] = 0.;
            }
            else
            {
                atm.ice_water_content[i] *= air_density[i]*g_per_kg;
            }
        }
        alloc(atm.liquid_water_content, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
        nc_catch(nc_inq_varid(ncid, "clwc", &varid));
        start[0] = t; start[1] = z; start[2] = y; start[3] = x;
        count[0] = atm.num_times; count[1] = atm.num_layers; count[2] = nlat; count[3] = nlon;
        get_var(ncid, varid, start, count, clouds);
        tzyx_to_tyxz(atm.liquid_water_content, clouds, nlon, nlat, atm.num_layers, atm.num_times);
        free(clouds);
        for (i=0; i<atm.num_times*atm.num_columns*atm.num_layers; ++i)
        {
            if (atm.liquid_water_content[i] <= 0.)
            {
                atm.liquid_water_content[i] = 0.;
            }
            else
            {
                atm.liquid_water_content[i] *= air_density[i]*g_per_kg;
            }
        }
        free(air_density);

        /*Calculate layer thickness.*/
        alloc(atm.layer_thickness, atm.num_times*atm.num_columns*atm.num_layers, fp_t *);
        for (i=0; i<atm.num_times; ++i)
        {
            int j;
            for (j=0; j<nlat; ++j)
            {
                int k;
                for (k=0; k<nlon; ++k)
                {
                    int offset = atm.num_levels*(i*nlat*nlon + j*nlon + k);
                    fp_t const *plev = &(atm.level_pressure[offset]);
                    offset = atm.num_layers*(i*nlat*nlon + j*nlon + k);
                    fp_t const *tlay = &(atm.layer_temperature[offset]);
                    double const gravity = 9.81; /*[m s-2].*/
                    int m;
                    for (m=0; m<atm.num_layers; ++m)
                    {
                        atm.layer_thickness[offset+m] = (fabs(log(plev[m]) - log(plev[m+1]))*
                                                        tlay[m]*gas_constant)/
                                                        (molar_mass*kg_per_g*gravity);
                    }
                }
            }
        }
    }

    /*Close the era5 file.*/
    nc_catch(nc_close(ncid));

    /*Surface emissivity. CIRC cases assume this is 1.*/
    atm.emissivity_grid_size = 2;
    alloc(atm.emissivity_grid, atm.emissivity_grid_size, fp_t *);
    atm.emissivity_grid[0] = -1.;
    atm.emissivity_grid[1] = 0.;
    alloc(atm.surface_emissivity, atm.num_times*atm.num_columns*atm.emissivity_grid_size, fp_t *);
    for (i=0; i<atm.num_times*atm.num_columns*atm.emissivity_grid_size; ++i)
    {
        atm.surface_emissivity[i] = 1.;
    }

    /*Open the greenhouse gas file.*/
    get_argument(*parser, "ghg_file", buffer);
    nc_catch(nc_open(buffer, NC_NOWRITE, &ncid));
    get_argument(*parser, "-year", buffer);
    int year = atoi(buffer);

    /*Greenhouse gas abundances.*/
    struct MoleculeMeta ghg_molecules[3] = {{CH4, "-CH4", "ch4", 0.},
                                            {CO2, "-CO2", "co2", 0.},
                                            {N2O, "-N2O", "n2o", 0.}};
    alloc(abundance, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    for (i=0; i<3; ++i)
    {
        if (get_argument(*parser, ghg_molecules[i].flag, NULL))
        {
            atm.molecules[atm.num_molecules] = ghg_molecules[i].id;
            alloc(atm.ppmv[atm.num_molecules], atm.num_times*atm.num_columns*atm.num_levels, fp_t *);
            fp_t *ppmv = atm.ppmv[atm.num_molecules];
            nc_catch(nc_inq_varid(ncid, ghg_molecules[i].name, &varid));
            start[0] = year - 1; start[1] = 0; start[2] = 0; start[3] = 0;
            count[0] = 1; count[1] = 1; count[2] = 1; count[3] = 1;
            fp_t ab_in;
            get_var(ncid, varid, start, count, &ab_in);
            int j;
            for (j=0; j<atm.num_times*atm.num_levels*nlat*nlon; ++j)
            {
                abundance[j] = ab_in;
            }
            tzyx_to_tyxz(ppmv, abundance, nlon, nlat, atm.num_levels, atm.num_times);
            atm.num_molecules++;
        }
    }
    free(abundance);

    /*CFC abundances.*/
    atm.num_cfcs = 0;
    int const num_cfcs = 4;
    struct MoleculeMeta cfcs[num_cfcs] = {{CFC11, "-CFC-11", "f11", 0.},
                                          {CFC12, "-CFC-12", "f12", 0.},
                                          {HCFC22, "-HCFC-22", "f22", 0.},
                                          {CFC113, "-CFC-113", "f113", 0.}};
    alloc(abundance, atm.num_times*atm.num_levels*nlat*nlon, fp_t *);
    alloc(atm.cfc, num_cfcs, Cfc_t *);
    alloc(atm.cfc_ppmv, num_cfcs, fp_t **);
    for (i=0; i<num_cfcs; ++i)
    {
        if (get_argument(*parser, cfcs[i].flag, atm.cfc[atm.num_cfcs].path))
        {
            atm.cfc[atm.num_cfcs].id = cfcs[i].id;
            alloc(atm.cfc_ppmv[atm.num_cfcs], atm.num_times*atm.num_levels*atm.num_columns, fp_t *);
            nc_catch(nc_inq_varid(ncid, cfcs[i].name, &varid));
            start[0] = year - 1; start[1] = 0; start[2] = 0; start[3] = 0;
            count[0] = 1; count[1] = 1; count[2] = 1; count[3] = 1;
            fp_t ab_in;
            get_var(ncid, varid, start, count, &ab_in);
            fp_t *ppmv = atm.cfc_ppmv[atm.num_cfcs];
            int j;
            for (j=0; j<atm.num_times*atm.num_levels*nlat*nlon; ++j)
            {
                abundance[j] = ab_in;
            }
            tzyx_to_tyxz(ppmv, abundance, nlon, nlat, atm.num_levels, atm.num_times);
            atm.num_cfcs++;
        }
    }
    free(abundance);

    /*Collision-induced absorption (CIA) abundances.*/
    atm.num_cias = 0;
    atm.num_cia_species = 0;
    int const num_cias = 3;
    int const num_cia_species = 2;
    struct CiaMeta
    {
        int species1;
        int species2;
        char *flag;
    };
    struct CiaMeta cias[num_cias] = {{CIA_N2, CIA_N2, "-N2-N2"},
                                     {CIA_O2, CIA_N2, "-O2-N2"},
                                     {CIA_O2, CIA_O2, "-O2-O2"}};
    alloc(atm.cia, num_cias, Cia_t *);
    alloc(atm.cia_species, num_cia_species, int *);
    alloc(atm.cia_ppmv, num_cia_species, fp_t **);
    for (i=0; i<num_cias; ++i)
    {
        if (get_argument(*parser, cias[i].flag, atm.cia[atm.num_cias].path))
        {
            atm.cia[atm.num_cias].id[0] = cias[i].species1;
            atm.cia[atm.num_cias].id[1] = cias[i].species2;
            int j;
            for (j=0; j<2; ++j)
            {
                int k;
                for (k=0; k<atm.num_cia_species; ++k)
                {
                    if (atm.cia[atm.num_cias].id[j] == atm.cia_species[k])
                    {
                        break;
                    }
                }
                if (k >= atm.num_cia_species)
                {
                    atm.cia_species[atm.num_cia_species] = atm.cia[atm.num_cias].id[j];
                    alloc(atm.cia_ppmv[atm.num_cia_species], atm.num_times*atm.num_levels*atm.num_columns, fp_t *);
                    fp_t *ppmv = atm.cia_ppmv[atm.num_cia_species];
                    if (atm.cia_species[atm.num_cia_species] == CIA_N2)
                    {
                        for (k=0; k<atm.num_times*atm.num_columns*atm.num_levels; ++k)
                        {
                            ppmv[k] = 0.781*to_ppmv;
                        }
                    }
                    else if (atm.cia_species[atm.num_cia_species] == CIA_O2)
                    {
                        for (k=0; k<atm.num_times*atm.num_columns*atm.num_levels; ++k)
                        {
                            ppmv[k] = 0.21*to_ppmv;
                        }
                    }
                    atm.num_cia_species++;
                }
            }
            atm.num_cias++;
        }
    }
    nc_catch(nc_close(ncid));

    /*Aerosols.*/
    atm.clean = 1;
    return atm;
}


/*Free memory for atmosphere.*/
void destroy_atmosphere(Atmosphere_t * const atm)
{
    free(atm->layer_pressure);
    free(atm->level_pressure);
    free(atm->layer_temperature);
    free(atm->level_temperature);
    free(atm->solar_zenith_angle);
    free(atm->surface_temperature);
    free(atm->total_solar_irradiance);
    free(atm->albedo_grid);
    free(atm->surface_albedo);
    free(atm->emissivity_grid);
    free(atm->surface_emissivity);
    if (!atm->clean)
    {
        free(atm->aerosol_grid);
        free(atm->aerosol_optical_depth);
        free(atm->aerosol_single_scatter_albedo);
        free(atm->aerosol_asymmetry_factor);
    }
    if (!atm->clear)
    {
        free(atm->cloud_fraction);
        free(atm->ice_water_content);
        free(atm->liquid_water_content);
        free(atm->layer_thickness);
    }
    int i;
    for (i=0; i<atm->num_molecules; ++i)
    {
        free(atm->ppmv[i]);
    }
    free(atm->ppmv);
    free(atm->molecules);
    for (i=0; i<atm->num_cfcs; ++i)
    {
        free(atm->cfc_ppmv[i]);
    }
    if (atm->num_cfcs > 0)
    {
        free(atm->cfc_ppmv);
        free(atm->cfc);
    }
    for (i=0; i<atm->num_cia_species; ++i)
    {
        free(atm->cia_ppmv[i]);
    }
    if (atm->num_cias > 0)
    {
        free(atm->cia_ppmv);
        free(atm->cia);
        free(atm->cia_species);
    }
    return;
}


/*Add a variable to the output file.*/
static void add_flux_variable(Output_t * const o, /*Output object.*/
                              VarId_t const index, /*Variable index.*/
                              char const * const name, /*Variable name.*/
                              char const * const standard_name, /*Variable standard name.*/
                              char const * const units, /*Units.*/
                              fp_t const * const fill_value, /*Fill value.*/
                              int const num_dims /*Number of dimensions.*/
                             )
{
#ifdef SINGLE_PRESCISION
    nc_type type = NC_FLOAT;
#else
    nc_type type = NC_DOUBLE;
#endif
    int varid;
    int dimids[5] = {o->dimid[TIME], o->dimid[LEVEL], o->dimid[LAT], o->dimid[LON],
                     o->dimid[LW_WAVENUMBER]};
    nc_catch(nc_def_var(o->ncid, name, type, num_dims, dimids, &varid));
    nc_catch(nc_put_att_text(o->ncid, varid, "units", strlen(units), units));
    nc_catch(nc_put_att_text(o->ncid, varid, "standard_name", strlen(standard_name),
                             standard_name));
    if (fill_value != NULL)
    {
#ifdef SINGLE_PRESCISION
        nc_catch(nc_put_att_float(o->ncid, varid, "_FillValue", type, 1, fill_value));
#else
        nc_catch(nc_put_att_double(o->ncid, varid, "_FillValue", type, 1, fill_value));
#endif
    }
    o->varid[index] = varid;
    return;
}


/*Create an output file and write metadata.*/
void create_flux_file(Output_t **output, char const * const filepath,
                      Atmosphere_t const * const atm, SpectralGrid_t const * const lw_grid,
                      SpectralGrid_t const * const sw_grid)
{
    Output_t *file = (Output_t *)malloc(sizeof(*file));
    file->dimid = (int *)malloc(sizeof(*(file->dimid))*NUM_DIMS);
    file->varid = (int *)malloc(sizeof(*(file->varid))*NUM_VARS);
    nc_catch(nc_create(filepath, NC_NETCDF4, &(file->ncid)));
    nc_catch(nc_put_att_int(file->ncid, NC_GLOBAL, "lon_start", NC_INT, 1, &x));
    nc_catch(nc_put_att_int(file->ncid, NC_GLOBAL, "lon_stop", NC_INT, 1, &X));
    int lon_size = (int)global_num_lon;
    nc_catch(nc_put_att_int(file->ncid, NC_GLOBAL, "lon_global_size", NC_INT, 1, &lon_size));

    nc_catch(nc_def_dim(file->ncid, "time", atm->num_times, &(file->dimid[TIME])));
    int varid;
    nc_catch(nc_def_var(file->ncid, "time", NC_DOUBLE, 1, &(file->dimid[TIME]), &varid));
    nc_catch(nc_put_att_text(file->ncid, varid, "units", strlen(time_units), time_units));
    nc_catch(nc_put_att_text(file->ncid, varid, "axis", 1, "T"));
    size_t start[1] = {0};
    size_t count[1] = {atm->num_times};
    nc_catch(nc_put_vara_double(file->ncid, varid, start, count, times));

    nc_catch(nc_def_dim(file->ncid, "level", atm->num_levels, &(file->dimid[LEVEL])));

    nc_catch(nc_def_dim(file->ncid, "lat", nlat, &(file->dimid[LAT])));
    nc_catch(nc_def_var(file->ncid, "lat", NC_DOUBLE, 1, &(file->dimid[LAT]), &varid));
    nc_catch(nc_put_att_text(file->ncid, varid, "units", strlen(lat_units), lat_units));
    nc_catch(nc_put_att_text(file->ncid, varid, "axis", 1, "Y"));
    start[0] = 0; count[0] = nlat;
    nc_catch(nc_put_vara_double(file->ncid, varid, start, count, lats));

    nc_catch(nc_def_dim(file->ncid, "lon", nlon, &(file->dimid[LON])));
    nc_catch(nc_def_var(file->ncid, "lon", NC_DOUBLE, 1, &(file->dimid[LON]), &varid));
    nc_catch(nc_put_att_text(file->ncid, varid, "units", strlen(lon_units), lon_units));
    nc_catch(nc_put_att_text(file->ncid, varid, "axis", 1, "X"));
    start[0] = 0; count[0] = nlon;
    nc_catch(nc_put_vara_double(file->ncid, varid, start, count, lons));

    nc_catch(nc_def_dim(file->ncid, "layer", atm->num_layers, &(file->dimid[LAYER])));

    nc_catch(nc_def_dim(file->ncid, "wavenumber", lw_grid->n, &(file->dimid[LW_WAVENUMBER])));
    nc_catch(nc_def_var(file->ncid, "wavenumber", NC_DOUBLE, 1,
             &(file->dimid[LW_WAVENUMBER]), &varid));
    char *units = "cm-1";
    nc_catch(nc_put_att_text(file->ncid, varid, "units", strlen(units), units));
    start[0] = 0; count[0] = lw_grid->n;
    double grid[lw_grid->n];
    size_t i;
    for (i=0; i<lw_grid->n; ++i)
    {
        grid[i] = lw_grid->w0 + i*lw_grid->dw;
    }
    nc_catch(nc_put_vara_double(file->ncid, varid, start, count, grid));

    fp_t const fill = -1;
/*
    add_flux_variable(file, RLUAF, "rluaf", "upwelling_aerosol_free_longwave_flux_in_air",
                      "W m-2", &fill, 5);
    add_flux_variable(file, RLUCSAF, "rlucsaf",
                      "upwelling_clear_sky_aerosol_free_longwave_flux_in_air", "W m-2",
                      &fill, 5);
    add_flux_variable(file, RLDAF, "rldaf", "downwelling_aerosol_free_longwave_flux_in_air",
                      "W m-2", &fill, 5);
    add_flux_variable(file, RLDCSAF, "rldcsaf",
                      "downwelling_clear_sky_aerosol_free_longwave_flux_in_air",
                      "W m-2", &fill, 5);
*/
    add_flux_variable(file, PLEV, "p", "air_pressure", "mb", NULL, 4);
    add_flux_variable(file, TLEV, "t", "air_temperature", "K", NULL, 4);
    add_flux_variable(file, H2OVMR, "h2o_vmr", "water_vapor_vmr", "ppmv", NULL, 4);
    add_flux_variable(file, O3VMR, "o3_vmr", "ozone_vmr", "ppmv", NULL, 4);
    add_flux_variable(file, CH4VMR, "ch4_vmr", "methane_vmr", "ppmv", NULL, 4);
    add_flux_variable(file, CO2VMR, "co2_vmr", "carbon_dioxide_vmr", "ppmv", NULL, 4);
    add_flux_variable(file, N2OVMR, "n2o_vmr", "nitrous_oxide_vmr", "ppmv", NULL, 4);
#ifdef SINGLE_PRESCISION
    nc_type type = NC_FLOAT;
#else
    nc_type type = NC_DOUBLE;
#endif
    int dimids[5] = {file->dimid[TIME], file->dimid[LAT], file->dimid[LON]};
    nc_catch(nc_def_var(file->ncid, "ts", type, 3, dimids, &(file->varid[TS])));
    char *standard_name = "surface_temperature";
    nc_catch(nc_put_att_text(file->ncid, file->varid[TS], "standard_name", strlen(standard_name),
                             standard_name));
    units = "K";
    nc_catch(nc_put_att_text(file->ncid, file->varid[TS], "units", strlen(units), units));

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAYER];
    dimids[2] = file->dimid[LAT]; dimids[3] = file->dimid[LON];
    nc_catch(nc_def_var(file->ncid, "t_layer", type, 4, dimids, &(file->varid[TLAY])));
    standard_name = "air_layer_temperature";
    nc_catch(nc_put_att_text(file->ncid, file->varid[TLAY], "standard_name", strlen(standard_name),
                             standard_name));
    units = "K";
    nc_catch(nc_put_att_text(file->ncid, file->varid[TLAY], "units", strlen(units), units));
/*
    dimids[0] = file->dimid[TIME]; dimids[2] = file->dimid[LAYER];
    dimids[2] = file->dimid[LAT]; dimids[3] = file->dimid[LON];
    dimids[4] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "taucsaf", type, 5, dimids, &(file->varid[TAU])));
    standard_name = "clear_sky_aerosol_free_optical_depth";
    nc_catch(nc_put_att_text(file->ncid, file->varid[TAU], "standard_name", strlen(standard_name),
                             standard_name));
*/

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAT]; dimids[2] = file->dimid[LON];
    dimids[3] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "tau_z", type, 4, dimids, &(file->varid[TAUZ])));
    standard_name = "clear_sky_aerosol_free_integrated_optical_depth";
    nc_catch(nc_put_att_text(file->ncid, file->varid[TAUZ], "standard_name", strlen(standard_name),
                             standard_name));

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAT]; dimids[2] = file->dimid[LON];
    dimids[3] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "rldsaf", type, 4, dimids, &(file->varid[RLDSAF])));
    standard_name = "downwelling_surface_aerosol_free_longwave_flux_in_air";
    nc_catch(nc_put_att_text(file->ncid, file->varid[RLDSAF], "standard_name",
                             strlen(standard_name), standard_name));

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAT]; dimids[2] = file->dimid[LON];
    dimids[3] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "rlusaf", type, 4, dimids, &(file->varid[RLUSAF])));
    standard_name = "upwelling_surface_aerosol_free_longwave_flux_in_air";
    nc_catch(nc_put_att_text(file->ncid, file->varid[RLUSAF], "standard_name",
                             strlen(standard_name), standard_name));

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAT]; dimids[2] = file->dimid[LON];
    dimids[3] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "rlutaf", type, 4, dimids, &(file->varid[RLUTAF])));
    standard_name = "upwelling_toa_aerosol_free_longwave_flux_in_air";
    nc_catch(nc_put_att_text(file->ncid, file->varid[RLUTAF], "standard_name",
                             strlen(standard_name), standard_name));

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAT]; dimids[2] = file->dimid[LON];
    dimids[3] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "rldscsaf", type, 4, dimids, &(file->varid[RLDSCSAF])));
    standard_name = "downwelling_surface_clear_sky_aerosol_free_longwave_flux_in_air";
    nc_catch(nc_put_att_text(file->ncid, file->varid[RLDSCSAF], "standard_name",
                             strlen(standard_name), standard_name));

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAT]; dimids[2] = file->dimid[LON];
    dimids[3] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "rluscsaf", type, 4, dimids, &(file->varid[RLUSCSAF])));
    standard_name = "upwelling_surface_clear_sky_aerosol_free_longwave_flux_in_air";
    nc_catch(nc_put_att_text(file->ncid, file->varid[RLUSCSAF], "standard_name",
                             strlen(standard_name), standard_name));

    dimids[0] = file->dimid[TIME]; dimids[1] = file->dimid[LAT]; dimids[2] = file->dimid[LON];
    dimids[3] = file->dimid[LW_WAVENUMBER];
    nc_catch(nc_def_var(file->ncid, "rlutcsaf", type, 4, dimids, &(file->varid[RLUTCSAF])));
    standard_name = "upwelling_toa_clear_sky_aerosol_free_longwave_flux_in_air";
    nc_catch(nc_put_att_text(file->ncid, file->varid[RLUTCSAF], "standard_name",
                             strlen(standard_name), standard_name));

    *output = file;
    return;
}


/*Close output file.*/
void close_flux_file(Output_t * const o)
{
    nc_catch(nc_close(o->ncid));
    free(o->varid);
    free(o->dimid);
    return;
}


/*Write fluxes to the output file.*/
void write_output(Output_t *output, VarId_t id, fp_t *data, int time, int column)
{
    int lat = column/nlon;
    int lon = column - nlon*lat;
    if (id == TS)
    {
        size_t start[3] = {time, lat, lon};
        size_t count[3] = {1, 1, 1};
#ifdef SINGLE_PRECISION
        nc_catch(nc_put_vara_float(output->ncid, output->varid[id], start, count, data));
#else
        nc_catch(nc_put_vara_double(output->ncid, output->varid[id], start, count, data));
#endif
    }
    else if (id == TLAY)
    {
        size_t num_layers;
        nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LAYER], &num_layers));
        size_t start[4] = {time, 0, lat, lon};
        size_t count[4] = {1, num_layers, 1, 1};
#ifdef SINGLE_PRECISION
        nc_catch(nc_put_vara_float(output->ncid, output->varid[id], start, count, data));
#else
        nc_catch(nc_put_vara_double(output->ncid, output->varid[id], start, count, data));
#endif
    }
    else if (id >= RLDAF && id <= RLUCSAF)
    {
        size_t num_levels;
        nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LEVEL], &num_levels));
        size_t num_w;
        nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LW_WAVENUMBER], &num_w));
        size_t start[5] = {time, 0, lat, lon, 0};
        size_t count[5] = {1, num_levels, 1, 1, num_w};
#ifdef SINGLE_PRECISION
        nc_catch(nc_put_vara_float(output->ncid, output->varid[id], start, count, data));
#else
        nc_catch(nc_put_vara_double(output->ncid, output->varid[id], start, count, data));
#endif
    }
    else if (id == TAU)
    {
        size_t num_layers;
        nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LAYER], &num_layers));
        size_t num_w;
        nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LW_WAVENUMBER], &num_w));
        size_t start[5] = {time, 0, lat, lon, 0};
        size_t count[5] = {1, num_layers, 1, 1, num_w};
#ifdef SINGLE_PRECISION
        nc_catch(nc_put_vara_float(output->ncid, output->varid[id], start, count, data));
#else
        nc_catch(nc_put_vara_double(output->ncid, output->varid[id], start, count, data));
#endif
    }
    else if (id == TAUZ || id == RLDSAF || id == RLUSAF || id == RLUTAF ||
             id == RLDSCSAF || id == RLUSCSAF || id == RLUTCSAF)
    {
        size_t num_w;
        nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LW_WAVENUMBER], &num_w));
        size_t start[4] = {time, lat, lon, 0};
        size_t count[4] = {1, 1, 1, num_w};
#ifdef SINGLE_PRECISION
        nc_catch(nc_put_vara_float(output->ncid, output->varid[id], start, count, data));
#else
        nc_catch(nc_put_vara_double(output->ncid, output->varid[id], start, count, data));
#endif
    }
    else
    {
        size_t num_levels;
        nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LEVEL], &num_levels));
        size_t start[4] = {time, 0, lat, lon};
        size_t count[4] = {1, num_levels, 1, 1};
#ifdef SINGLE_PRECISION
        nc_catch(nc_put_vara_float(output->ncid, output->varid[id], start, count, data));
#else
        nc_catch(nc_put_vara_double(output->ncid, output->varid[id], start, count, data));
#endif
    }
    return;
}
