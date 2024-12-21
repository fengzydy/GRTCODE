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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "argparse.h"
#include "clouds_lib.h"
#include "driver.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"
#include "longwave.h"
#include "rayleigh.h"
#include "shortwave.h"
#include "solar_flux.h"


#define catch(e) { \
    int e_ = e; \
    if (e_ != GRTCODE_SUCCESS) { \
        fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
        char b_[1024]; \
        grtcode_errstr(e_, b_, 1024); \
        fprintf(stderr, "%s", b_); \
        return EXIT_FAILURE; \
    }}

#define LONGWAVE_PASS 1
#define MAX_SIZE 64
#define SHORTWAVE_PASS 2

/*Atmospheric column.*/
typedef struct AtmosphericColumn {
    int clean;
    int clear;
    int num_cfcs;
    int num_cia_species;
    int num_layers;
    int num_levels;
    int num_molecules;
    int molecules[MAX_SIZE];
    int cfcs[MAX_SIZE];
    int cia_species[MAX_SIZE];
    fp_t cosine_diffuse_angle;
    fp_t cosine_zenith_angle;
    fp_t * diffuse_albedo;
    fp_t * direct_albedo;
    fp_t * emissivity;
    fp_t * layer_temperature;
    fp_t * level_pressure;
    fp_t * level_temperature;
    fp_t * overlap_parameter;
    fp_t surface_temperature;
    fp_t total_solar_irradiance;
    fp_t * cloud_fraction;
    fp_t * liquid_content;
    fp_t * ice_content;
    fp_t * thickness;
    fp_t * incident_solar_flux;
    fp_t const * ppmv[MAX_SIZE];
    fp_t const * cfc_ppmv[MAX_SIZE];
    fp_t const * cia_ppmv[MAX_SIZE];
} AtmosphericColumn_t;


/*Extract the data for an atmospheric column.*/
static int atmospheric_column(AtmosphericColumn_t * column,
                              SpectralGrid_t lw_grid,
                              SpectralGrid_t sw_grid,
                              Atmosphere_t atmosphere,
                              fp_t * surface_emissivity,
                              fp_t * direct_albedo,
                              fp_t * overlap,
                              fp_t * incident_solar_flux,
                              int offset,
                              int column_index
)
{
    int i = offset + column_index;
    column->num_layers = atmosphere.num_layers;
    column->num_levels = atmosphere.num_levels;
    column->clean = atmosphere.clean;
    column->clear = atmosphere.clear;

    /*Longwave radiation.*/
    int n = i*atmosphere.emissivity_grid_size;
    catch(interpolate_to_grid(lw_grid, atmosphere.emissivity_grid,
                              &(atmosphere.surface_emissivity[n]),
                              atmosphere.emissivity_grid_size, surface_emissivity,
                              linear_sample, constant_extrapolation));
    column->emissivity = surface_emissivity;
    column->surface_temperature = atmosphere.surface_temperature[i];

    /*Shortwave radiation.*/
    column->cosine_diffuse_angle = 0.5;
    column->cosine_zenith_angle = atmosphere.solar_zenith_angle[i];
    n = i*atmosphere.albedo_grid_size;
    catch(interpolate_to_grid(sw_grid, atmosphere.albedo_grid, &(atmosphere.surface_albedo[n]),
                              atmosphere.albedo_grid_size, direct_albedo,
                              linear_sample, constant_extrapolation));
    column->direct_albedo = direct_albedo;
    column->diffuse_albedo = direct_albedo;
    column->total_solar_irradiance = atmosphere.total_solar_irradiance[i];
    column->layer_temperature = &(atmosphere.layer_temperature[i*atmosphere.num_layers]);
    column->level_temperature = &(atmosphere.level_temperature[i*atmosphere.num_levels]);
    column->incident_solar_flux = incident_solar_flux;

    /*Shared.*/
    column->level_pressure = &(atmosphere.level_pressure[i*atmosphere.num_levels]);
    column->num_molecules = atmosphere.num_molecules;
    if (atmosphere.num_molecules > MAX_SIZE)
    {
        fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
        fprintf(stderr, "MAX_SIZE is too small for the number of molecules (%d) used.",
                atmosphere.num_molecules);
            return EXIT_FAILURE;
    }
    int j;
    for (j=0; j<atmosphere.num_molecules; ++j)
    {
        column->molecules[j] = atmosphere.molecules[j];
        fp_t const * ppmv = atmosphere.ppmv[j];
        column->ppmv[j] = &(ppmv[i*atmosphere.num_levels]);
    }
    if (atmosphere.num_cfcs > MAX_SIZE)
    {
        fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
        fprintf(stderr, "MAX_SIZE is too small for the number of cfcs (%d) used.",
                atmosphere.num_cfcs);
            return EXIT_FAILURE;
    }
    column->num_cfcs = atmosphere.num_cfcs;
    for (j=0; j<atmosphere.num_cfcs; ++j)
    {
        column->cfcs[j] = atmosphere.cfc[j].id;
        fp_t const * ppmv = atmosphere.cfc_ppmv[j];
        column->cfc_ppmv[j] = &(ppmv[i*atmosphere.num_levels]);
    }
    if (atmosphere.num_cfcs > MAX_SIZE)
    {
        fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
        fprintf(stderr, "MAX_SIZE is too small for the number of CIA species (%d) used.",
                atmosphere.num_cia_species);
            return EXIT_FAILURE;
    }
    column->num_cia_species = atmosphere.num_cia_species;
    for (j=0; j<atmosphere.num_cia_species; ++j)
    {
        column->cia_species[j] = atmosphere.cia_species[j];
        fp_t const * ppmv = atmosphere.cia_ppmv[j];
        column->cia_ppmv[j] = &(ppmv[i*atmosphere.num_levels]);
    }
    if (!atmosphere.clear)
    {
        /*Calculate overlap for the column.*/
        fp_t const * pressure = &(atmosphere.layer_pressure[i*atmosphere.num_layers]);
        fp_t altitude[atmosphere.num_layers];
        fp_t const pa_per_mb = 100.;
        fp_t const pressure_scale_height = 7.3; /*[km].*/
        int j;
        for (j=0; j<atmosphere.num_layers; ++j)
        {
            altitude[j] = log(pa_per_mb*pressure[j])*pressure_scale_height;
        }
        fp_t const scale_length = 2.;
        calculate_overlap(atmosphere.num_layers, altitude, scale_length, overlap);
        column->overlap_parameter = overlap;
        column->cloud_fraction = &(atmosphere.cloud_fraction[i*atmosphere.num_layers]);
        column->liquid_content = &(atmosphere.liquid_water_content[i*atmosphere.num_layers]);
        column->ice_content = &(atmosphere.ice_water_content[i*atmosphere.num_layers]);
        column->thickness = &(atmosphere.layer_thickness[i*atmosphere.num_layers]);
    }
    return GRTCODE_SUCCESS;
}


/*Add molecules, CFCs, and collision-induced absoprtion.*/
static int add_molecules(GasOptics_t * const lbl, /*Gas optics object.*/
                         Atmosphere_t const atm /*Atmospheric state.*/
)
{
    int i;
    for (i=0; i<atm.num_molecules; ++i)
    {
        catch(add_molecule(lbl, atm.molecules[i], NULL, NULL));
    }
    for (i=0; i<atm.num_cfcs; ++i)
    {
        catch(add_cfc(lbl, atm.cfc[i].id, atm.cfc[i].path));
    }
    for (i=0; i<atm.num_cias; ++i)
    {
        catch(add_cia(lbl, atm.cia[i].id[0], atm.cia[i].id[1], atm.cia[i].path));
    }
    return GRTCODE_SUCCESS;
}


/*Calculate aerosol optics.*/
static int calculate_aerosol_optics(AtmosphericColumn_t atm_column, /*Atmospheric state.*/
                                    Optics_t * const aerosol /*Optics due to aerosols.*/
)
{
    int n = atm_column.num_layers*aerosol->grid.n;
    fp_t * tau = (fp_t *)malloc(sizeof(*tau)*n);
    fp_t * omega = (fp_t *)malloc(sizeof(*omega)*n);
    fp_t * g = (fp_t *)malloc(sizeof(*g)*n);
/*
    int j;
    for (j=0; j<atm_column.num_layers; ++j)
    {
        n = i*atm_column.num_layers*atm_column.aerosol_grid_size + j*atm_column.aerosol_grid_size;
        catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_optical_depth[n]),
                                  atm.aerosol_grid_size, &(tau[j*aerosol->grid.n]),
                                  linear_sample, NULL));
        catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_single_scatter_albedo[n]),
                                  atm.aerosol_grid_size, &(omega[j*aerosol->grid.n]),
                                  linear_sample, NULL));
        catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_asymmetry_factor[n]),
                                  atm.aerosol_grid_size, &(g[j*aerosol->grid.n]),
                                  linear_sample, NULL));
    }
    catch(update_optics(aerosol, tau, omega, g));
*/
    free(tau);
    free(omega);
    free(g);
    return GRTCODE_SUCCESS;
}


/*Calculate gas optics.*/
static int calculate_gas_optics(GasOptics_t * const lbl, /*Gas optics object.*/
                                AtmosphericColumn_t const atm_column, /*Atmospheric state.*/
                                Optics_t * const gas, /*Optics due to molecular lines.*/
                                Optics_t * const rayleigh /*Optics due to rayleigh scattering.*/
)
{
    int i;
    for (i=0; i<atm_column.num_molecules; ++i)
    {
        catch(set_molecule_ppmv(lbl, atm_column.molecules[i], atm_column.ppmv[i]));
    }
    for (i=0; i<atm_column.num_cfcs; ++i)
    {
        catch(set_cfc_ppmv(lbl, atm_column.cfcs[i], atm_column.cfc_ppmv[i]));
    }
    for (i=0; i<atm_column.num_cia_species; ++i)
    {
        catch(set_cia_ppmv(lbl, atm_column.cia_species[i], atm_column.cia_ppmv[i]));
    }
    catch(calculate_optical_depth(lbl, atm_column.level_pressure,
                                  atm_column.level_temperature, gas));
    catch(rayleigh_scattering(rayleigh, atm_column.level_pressure));
    return GRTCODE_SUCCESS;
}


/*Calculate the fluxes in a single atmospheric column.*/
static int column_calculation(int label,
                              AtmosphericColumn_t atm_column,
                              GasOptics_t lbl,
                              Optics_t optics_gas,
                              Optics_t optics_rayleigh,
                              Optics_t optics_aerosol,
                              Optics_t optics_liquid_cloud,
                              Optics_t optics_ice_cloud,
                              void * solver_object,
                              fp_t * flux_up,
                              fp_t * flux_down,
                              SpectralGrid_t grid,
                              Output_t * output,
                              int time_index,
                              int column_index
)
{
    /*Gas optics.*/
    Optics_t optics_total;
    catch(calculate_gas_optics(&lbl, atm_column, &optics_gas, &optics_rayleigh));
    Optics_t const * optics_array[4] = {&optics_gas, &optics_rayleigh, NULL, NULL};
    catch(add_optics(optics_array, 2, &optics_total));

    unsigned int regime;
    if (label == LONGWAVE_PASS)
    {
        /*Longwave clear-clean-sky fluxes.*/
        catch(calculate_lw_fluxes((Longwave_t *)solver_object, &optics_total,
                                  atm_column.surface_temperature,
                                  atm_column.layer_temperature,
                                  atm_column.level_temperature,
                                  atm_column.emissivity,
                                  flux_up, flux_down));
        regime = LONGWAVE;
    }
    else
    {
        /*Shortwave clear-clean-sky fluxes.*/
        catch(calculate_sw_fluxes((Shortwave_t *)solver_object, &optics_total,
                                  atm_column.cosine_zenith_angle,
                                  atm_column.cosine_diffuse_angle,
                                  atm_column.direct_albedo,
                                  atm_column.diffuse_albedo,
                                  atm_column.total_solar_irradiance,
                                  atm_column.incident_solar_flux,
                                  flux_up, flux_down));
        regime = SHORTWAVE;
    }
    catch(destroy_optics(&optics_total));
    int n = grid.n*(atm_column.num_levels - 1);
    write_output(output, (unsigned int)(regime ^ UPWARD ^ TOP ^ CLEARSKY ^ AEROSOLFREE),
                 flux_up, time_index, column_index);
    write_output(output, (unsigned int)(regime ^ UPWARD ^ SURFACE ^ CLEARSKY ^ AEROSOLFREE),
                 &(flux_up[n]), time_index, column_index);
    write_output(output, (unsigned int)(regime ^ DOWNWARD ^ SURFACE ^ CLEARSKY ^ AEROSOLFREE),
                 &(flux_down[n]), time_index, column_index);
    if (label != LONGWAVE_PASS)
    {
        write_output(output, (unsigned int)(regime ^ DOWNWARD ^ TOP ^ CLEARSKY ^ AEROSOLFREE),
                     flux_down, time_index, column_index);
    }

    if (!atm_column.clean)
    {
        /*Aerosol optics.*/
        catch(calculate_aerosol_optics(atm_column, &optics_aerosol));
        optics_array[2] = &optics_aerosol;
        catch(add_optics(optics_array, 3, &optics_total));

        if (label == LONGWAVE_PASS)
        {
            /*Longwave clear-sky fluxes.*/
            catch(calculate_lw_fluxes((Longwave_t *)solver_object, &optics_total,
                                      atm_column.surface_temperature,
                                      atm_column.layer_temperature,
                                      atm_column.level_temperature,
                                      atm_column.emissivity,
                                      flux_up, flux_down));
        }
        else
        {
            /*Shortwave clear-sky fluxes.*/
            catch(calculate_sw_fluxes((Shortwave_t *)solver_object, &optics_total,
                                      atm_column.cosine_zenith_angle,
                                      atm_column.cosine_diffuse_angle,
                                      atm_column.direct_albedo,
                                      atm_column.diffuse_albedo,
                                      atm_column.total_solar_irradiance,
                                      atm_column.incident_solar_flux,
                                      flux_up, flux_down));
        }
        catch(destroy_optics(&optics_total));
        n = grid.n*(atm_column.num_levels - 1);
        write_output(output, (unsigned int)(regime ^ UPWARD ^ TOP ^ CLEARSKY ^ AEROSOL),
                     flux_up, time_index, column_index);
        write_output(output, (unsigned int)(regime ^ UPWARD ^ SURFACE ^ CLEARSKY ^ AEROSOL),
                     &(flux_up[n]), time_index, column_index);
        write_output(output, (unsigned int)(regime ^ DOWNWARD ^ SURFACE ^ CLEARSKY ^ AEROSOL),
                     &(flux_down[n]), time_index, column_index);
    }

    if (!atm_column.clear)
    {
        /*Set up band structure arrays for the clouds library.*/
        fp_t band_centers[grid.n];
        uint64_t j;
        for (j=0; j<grid.n; ++j)
        {
            band_centers[j] = grid.w0 + j*grid.dw;
        }
        fp_t band_limits[grid.n + 1];
        for (j=1; j<grid.n; ++j)
        {
            band_limits[j] = 0.5*(band_centers[j - 1] + band_centers[j]);
        }
        band_limits[0] = band_centers[0] - grid.dw;
        if (band_limits[0] < 0.)
        {
            band_limits[0] = 0;
        }
        band_limits[grid.n] = band_centers[grid.n - 1] + grid.dw;

        /*Zero out the flux summing variables.*/
        fp_t flux_up_sum[atm_column.num_levels*grid.n];
        fp_t flux_down_sum[atm_column.num_levels*grid.n];
        for (j=0; j<atm_column.num_levels*grid.n; ++j)
        {
            flux_up_sum[j] = 0.;
            flux_down_sum[j] = 0.;
        }

        uint64_t const num_subcolumns = 5;
        for (j=0; j<num_subcolumns; ++j)
        {
            /*Cloud optics.*/
            cloud_optics((int)grid.n, band_centers, band_limits, atm_column.num_layers,
                         atm_column.cloud_fraction, atm_column.liquid_content,
                         atm_column.ice_content, atm_column.overlap_parameter, 10.,
                         atm_column.layer_temperature, optics_liquid_cloud.tau, 
                         optics_liquid_cloud.omega, optics_liquid_cloud.g,
                         optics_ice_cloud.tau, optics_ice_cloud.omega,
                         optics_ice_cloud.g);

            /*Convert from extinction coefficient to optical depth.*/
            uint64_t k;
            for (k=0; (int)k<atm_column.num_layers; ++k)
            {
                uint64_t m;
                for (m=0; m<grid.n; ++m)
                {
                    optics_liquid_cloud.tau[k*grid.n + m] *= atm_column.thickness[k];
                    optics_ice_cloud.tau[k*grid.n + m] *= atm_column.thickness[k];
                }
            }

            /*Add to gas optics.*/
            optics_array[2] = &optics_liquid_cloud;
            optics_array[3] = &optics_ice_cloud;
            catch(add_optics(optics_array, 4, &optics_total));

            if (label == LONGWAVE_PASS)
            {
                /*Longwave aerosol-free fluxes.*/
                catch(calculate_lw_fluxes((Longwave_t *)solver_object, &optics_total,
                                          atm_column.surface_temperature,
                                          atm_column.layer_temperature,
                                          atm_column.level_temperature,
                                          atm_column.emissivity,
                                          flux_up, flux_down));
            }
            else
            {
                /*Shortwave aerosol-free fluxes.*/
                catch(calculate_sw_fluxes((Shortwave_t *)solver_object, &optics_total,
                                          atm_column.cosine_zenith_angle,
                                          atm_column.cosine_diffuse_angle,
                                          atm_column.direct_albedo,
                                          atm_column.diffuse_albedo,
                                          atm_column.total_solar_irradiance,
                                          atm_column.incident_solar_flux,
                                          flux_up, flux_down));
            }
            catch(destroy_optics(&optics_total));
            for (k=0; k<atm_column.num_levels*grid.n; ++k)
            {
                flux_up_sum[k] += flux_up[k];
                flux_down_sum[k] += flux_down[k];
            }
        }
        for (j=0; j<atm_column.num_levels*grid.n; ++j)
        {
            flux_up_sum[j] /= (double)num_subcolumns;
            flux_down_sum[j] /= (double)num_subcolumns;
        }
        n = grid.n*(atm_column.num_levels - 1);
        write_output(output, (unsigned int)(regime ^ UPWARD ^ TOP ^ CLOUDYSKY ^ AEROSOLFREE),
                     flux_up, time_index, column_index);
        write_output(output, (unsigned int)(regime ^ UPWARD ^ SURFACE ^ CLOUDYSKY ^ AEROSOLFREE), 
                     &(flux_up[n]), time_index, column_index);
        write_output(output, (unsigned int)(regime ^ DOWNWARD ^ SURFACE ^ CLOUDYSKY ^ AEROSOLFREE),
                     &(flux_down[n]), time_index, column_index);
    }
    return GRTCODE_SUCCESS;
}


static int driver(Atmosphere_t const atm, /*Atmospheric state.*/
                  char const * const hitran_path, /*Path to HITRAN database ascii file.*/
                  char const * const solar_flux_path, /*Path to solar flux ascii file.*/
                  SpectralGrid_t const lw_grid, /*Longwave spectral grid.*/
                  SpectralGrid_t const sw_grid, /*Shortwave spectral grid.*/
                  Device_t const device, /*Device object.*/
                  Output_t * const output, /*Output object.*/
                  char const * const beta_path, /*Path to beta distribution input file.*/
                  char const * const ice_path, /*Path to ice cloud parameterization input file.*/
                  char const * const liquid_path /*Path to liquid cloud parameterization input file.*/
                 )
{
    /*Intialize gas optics objects.*/
    GasOptics_t lbl_lw;
    int method = line_sample;
    catch(create_gas_optics(&lbl_lw, atm.num_levels, &lw_grid, &device,
                            hitran_path, atm.h2o_ctm, atm.o3_ctm, NULL, &method));
    catch(add_molecules(&lbl_lw, atm));
    GasOptics_t lbl_sw;
    catch(create_gas_optics(&lbl_sw, atm.num_levels, &sw_grid, &device,
                            hitran_path, atm.h2o_ctm, atm.o3_ctm, NULL, &method));
    catch(add_molecules(&lbl_sw, atm));

    /*Initialize optics objects.*/
    /*Gas.*/
    Optics_t optics_lw_gas;
    catch(create_optics(&optics_lw_gas, atm.num_layers, &lw_grid, &device));
    Optics_t optics_lw_rayleigh;
    catch(create_optics(&optics_lw_rayleigh, atm.num_layers, &lw_grid, &device));
    Optics_t optics_sw_gas;
    catch(create_optics(&optics_sw_gas, atm.num_layers, &sw_grid, &device));
    Optics_t optics_sw_rayleigh;
    catch(create_optics(&optics_sw_rayleigh, atm.num_layers, &sw_grid, &device));

    /*Aerosol.*/
    Optics_t optics_lw_aerosol;
    Optics_t optics_sw_aerosol;
    if (!atm.clean)
    {
        catch(create_optics(&optics_lw_aerosol, atm.num_layers, &lw_grid, &device));
        catch(create_optics(&optics_sw_aerosol, atm.num_layers, &sw_grid, &device));
    }

    /*Clouds.*/
    Optics_t optics_lw_liquid_cloud;
    Optics_t optics_lw_ice_cloud;
    Optics_t optics_sw_liquid_cloud;
    Optics_t optics_sw_ice_cloud;
    if (!atm.clear)
    {
        catch(create_optics(&optics_lw_liquid_cloud, atm.num_layers, &lw_grid, &device));
        catch(create_optics(&optics_lw_ice_cloud, atm.num_layers, &lw_grid, &device));
        catch(create_optics(&optics_sw_liquid_cloud, atm.num_layers, &sw_grid, &device));
        catch(create_optics(&optics_sw_ice_cloud, atm.num_layers, &sw_grid, &device));

        /*Initialize clouds library.*/
        if (beta_path[0] == '\0' || ice_path[0] == '\0' || liquid_path[0] == '\0')
        {
            fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
            fprintf(stderr, "-beta-path, -ice-path, and -liquid-path args required when"
                            " running with clouds.");
            return EXIT_FAILURE;
        }
        initialize_clouds_lib(beta_path, ice_path, liquid_path, NULL);
    }

    /*Initialize the solar fluxe object.*/
    SolarFlux_t solar_flux;
    catch(create_solar_flux(&solar_flux, &sw_grid, solar_flux_path));

    /*Initialize solver objects.*/
    Longwave_t longwave;
    catch(create_longwave(&longwave, atm.num_levels, &lw_grid, &device));
    Shortwave_t shortwave;
    catch(create_shortwave(&shortwave, atm.num_levels, &sw_grid, &device));

    /*Initialize other buffers.*/
    fp_t * surface_emissivity = (fp_t *)malloc(sizeof(*surface_emissivity)*lw_grid.n);
    fp_t * lw_flux_up = (fp_t *)malloc(sizeof(*lw_flux_up)*atm.num_levels*lw_grid.n);
    fp_t * lw_flux_down = (fp_t *)malloc(sizeof(*lw_flux_down)*atm.num_levels*lw_grid.n);
    fp_t * direct_albedo = (fp_t *)malloc(sizeof(*direct_albedo)*sw_grid.n);
    fp_t * sw_flux_up = (fp_t *)malloc(sizeof(*sw_flux_up)*atm.num_levels*sw_grid.n);
    fp_t * sw_flux_down = (fp_t *)malloc(sizeof(*sw_flux_down)*atm.num_levels*sw_grid.n);
    fp_t * overlap = (fp_t *)malloc(sizeof(*overlap)*(atm.num_layers - 1));

    /*Loop through the columns.*/
    int t;
    for (t=0; t<atm.num_times; ++t)
    {
        int offset = t*atm.num_columns;
        int i;
        for (i=0; i<atm.num_columns; ++i)
        {
            AtmosphericColumn_t atm_column;
            catch(atmospheric_column(&atm_column, lw_grid, sw_grid, atm, surface_emissivity,
                                     direct_albedo, overlap, solar_flux.incident_flux,
                                     offset, i));
            catch(column_calculation(LONGWAVE_PASS, atm_column, lbl_lw, optics_lw_gas,
                                     optics_lw_rayleigh, optics_lw_aerosol,
                                     optics_lw_liquid_cloud, optics_lw_ice_cloud,
                                     (void *)&longwave, lw_flux_up, lw_flux_down, lw_grid,
                                     output, t, i));
            if (atm_column.cosine_zenith_angle > 0.)
            {
                catch(column_calculation(SHORTWAVE_PASS, atm_column, lbl_sw, optics_sw_gas,
                                         optics_sw_rayleigh, optics_sw_aerosol,
                                         optics_sw_liquid_cloud, optics_sw_ice_cloud,
                                         (void *)&shortwave, sw_flux_up, sw_flux_down,
                                         sw_grid, output, t, i));
            }

/*
            write_output(output, PLEV, atm_column.level_pressure, t, i);
            write_output(output, TLEV, atm_column.level_temperature, t, i);
            write_output(output, TLAY, atm_column.layer_temperature, t, i);
            write_output(output, TS, &(atm_column.surface_temperature), t, i);
            int j;
            for (j=0; j<atm.num_molecules; ++j)
            {
                write_output(output, H2OVMR + j, atm_column.ppmv[j], t, i);
            }
*/
        }
    }

    /*Release memory.*/
    catch(destroy_gas_optics(&lbl_lw));
    catch(destroy_gas_optics(&lbl_sw));
    catch(destroy_optics(&optics_lw_gas));
    catch(destroy_optics(&optics_lw_rayleigh));
    catch(destroy_optics(&optics_sw_gas));
    catch(destroy_optics(&optics_sw_rayleigh));
    if (!atm.clean)
    {
        catch(destroy_optics(&optics_lw_aerosol));
        catch(destroy_optics(&optics_sw_aerosol));
    }
    if (!atm.clear)
    {
        finalize_clouds_lib();
        catch(destroy_optics(&optics_lw_liquid_cloud));
        catch(destroy_optics(&optics_lw_ice_cloud));
        catch(destroy_optics(&optics_sw_liquid_cloud));
        catch(destroy_optics(&optics_sw_ice_cloud));
        initialize_clouds_lib(beta_path, ice_path, liquid_path, NULL);
    }
    catch(destroy_solar_flux(&solar_flux));
    catch(destroy_longwave(&longwave));
    catch(destroy_shortwave(&shortwave));
    free(surface_emissivity);
    free(lw_flux_up);
    free(lw_flux_down);
    free(direct_albedo);
    free(sw_flux_up);
    free(sw_flux_down);
    free(overlap);
    return GRTCODE_SUCCESS;
}


/*Main driver program.  When linking an executable, the user must provide an object file that
  includes implementations for the following defined in driver.h:
   - struct output;
   - void close_flux_file(Output_t *output);
   - Atmosphere_t create_atmosphere(Parser_t parser);
   - void destroy_atmosphere(Atmosphere_t *atm);
   - Output_t create_flux_file(char *path, Atmosphere_t *atm);
   - void write_output(Output_t *output, enum varid id, fp_t *data, int time, int column);*/
int main(int argc, char **argv)
{
    /*Add/parse command line arguments.*/
    char *description = "Calculates radiative fluxes using the line-by-line method.";
    Parser_t parser = create_parser(argc, argv, description);
    add_argument(&parser, "hitran_file", NULL, "HITRAN database file.", NULL);
    add_argument(&parser, "solar_flux", NULL, "Solar flux file.", NULL);
    int one = 1;
    add_argument(&parser, "-beta-path", NULL, "Path to beta distribution input file.", &one);
    add_argument(&parser, "-c", "--line-cutoff", "Cutoff [1/cm] from line center.", &one);
    add_argument(&parser, "-d", "--device", "GPU id", &one);
    add_argument(&parser, "-ice-path", NULL, "Path to ice cloud parameterization input file.", &one);
    add_argument(&parser, "-liquid-path", NULL, "Path to liquid cloud parameterization input file.", &one);
    add_argument(&parser, "-o", NULL, "Name of output file.", &one);
    add_argument(&parser, "-r-lw", "--lw-resolution", "Longwave spectral resolution [1/cm].", &one);
    add_argument(&parser, "-r-sw", "--sw-resolution", "Shortwave spectral resolution [1/cm].", &one);
    add_argument(&parser, "-s", "--solver", "Shortwave solver.", &one);
    add_argument(&parser, "-v", "--verbose", "Increase verbosity.", NULL);
    add_argument(&parser, "-w-lw", "--lw-lower-bound", "Longwave spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-w-sw", "--sw-lower-bound", "Shortwave spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-W-lw", "--lw-upper-bound", "Longwave spectral upper bound [1/cm].", &one);
    add_argument(&parser, "-W-sw", "--sw-upper-bound", "Shortwave spectral upper bound [1/cm].", &one);

    /*Get the atmospheric data.*/
    Atmosphere_t atm = create_atmosphere(&parser);

    /*Set GRTCODE verbosity.*/
    int verbosity_level = get_argument(parser, "-v", NULL) ? GRTCODE_INFO : GRTCODE_WARN;
    grtcode_set_verbosity(verbosity_level);

    /*Get paths to required input files.*/
    char hitran_path[valuelen];
    get_argument(parser, "hitran_file", hitran_path);
    char solar_flux_path[valuelen];
    get_argument(parser, "solar_flux", solar_flux_path);

    /*Create a spectral grids.*/
    char buffer[valuelen];
    double w0 = get_argument(parser, "-w-lw", buffer) ? atof(buffer) : 1.;
    double wn = get_argument(parser, "-W-lw", buffer) ? atof(buffer) : 3250.;
    double dw = get_argument(parser, "-r-lw", buffer) ? atof(buffer) : 0.1;
    SpectralGrid_t lw_grid;
    catch(create_spectral_grid(&lw_grid, w0, wn, dw));
    dw = get_argument(parser, "-r-sw", buffer) ? atof(buffer) : 1.;
    w0 = get_argument(parser, "-w-sw", buffer) ? atof(buffer) : 1.;
    wn = get_argument(parser, "-W-sw", buffer) ? atof(buffer) : 50000.;
    SpectralGrid_t sw_grid;
    catch(create_spectral_grid(&sw_grid, w0, wn, dw));

    /*Set the device to run on.*/
    Device_t device;
    int * device_id = NULL;
    if (get_argument(parser, "-d", buffer))
    {
        int d = atoi(buffer);
        device_id = &d;
    }
    catch(create_device(&device, device_id));

    /*Initialize the output file.*/
    if (!get_argument(parser, "-o", buffer))
    {
        snprintf(buffer, valuelen, "%s", "output.nc");
    }
    Output_t * output;
    create_flux_file(&output, buffer, &atm, &lw_grid, &sw_grid);

    /*Get cloud parameterization inputs.*/
    char beta_path[valuelen];
    if (!get_argument(parser, "-beta-path", beta_path))
    {
        beta_path[0] = '\0';
    }
    char ice_path[valuelen];
    if (!get_argument(parser, "-ice-path", ice_path))
    {
        ice_path[0] = '\0';
    }
    char liquid_path[valuelen];
    if (!get_argument(parser, "-liquid-path", liquid_path))
    {
        liquid_path[0] = '\0';
    }

    /*Calculate the fluxes over all the columns.*/
    catch(driver(atm, hitran_path, solar_flux_path, lw_grid, sw_grid, device, output,
                 beta_path, ice_path, liquid_path));

    /*Clean up.*/
    close_flux_file(output);
    destroy_atmosphere(&atm);
    destroy_parser(&parser);
    return EXIT_SUCCESS;
}
