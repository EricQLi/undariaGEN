// -------------------------------------------------------------------------- //
// organism.cpp                                                               //
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// UndariaGEN - A spatially-explicit agent-based model of invasive algal 
//            species in the marine environment
// Copyright (C) 2005-2015  James Thomas Murphy, Ray Walshe
//
// This file is a part of UndariaGEN
//
// UndariaGEN is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// UndariaGEN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
// Libraries                                                                  //
// -------------------------------------------------------------------------- //

#include "organism.hpp"
#include <stdlib.h>
#include <math.h>
#include<iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "bait.hpp"
#include <gsl/gsl_randist.h>


// -------------------------------------------------------------------------- //
// Static member initialisation                                               //
// -------------------------------------------------------------------------- //

const BaitParameters* Organism::baitParamsPtr_ = 0 ;
const EngineParameters* Organism::engParamsPtr_ = 0 ;
World<Patch> * Organism::ptrWorld = 0 ;
Fabric<Organism>* Organism::fabricPtr_         = 0 ;
AgentType_t Organism::agentType_ = ORGANISM;

int Organism::numSeedlings = 0;
int Organism::sumAgeMaturity = 0;
int Organism::sumAgeDeath = 0;
int Organism::countMature = 0;
int Organism::countDead = 0;
float Organism::gametGrowth=0.0;
float Organism::sporoGrowthMod=0.0;
double Organism::probGamMature=0.0;

float Organism::sumSporoGrowthMod=0.0;
unsigned Organism::countSporoGrowthMod=0;
float Organism::avgGametGrowthMod=0.0;
float Organism::avgSporoGrowthMod=0.0;


// -------------------------------------------------------------------------- //
// Constructor                                                                //
// -------------------------------------------------------------------------- //

Organism::Organism()
: Agent()
{
   stock_ = 0.0 ;
   age = 0;
   ageRecruit = 0;
   ageGameto = 0;
   ageMaturity = 0;
   sporeStock = 0.0;
   species=U_PINNATIFIDA;
   gametophyte=false;
   gamMature=false;
   tempType=NONE_T;
   shuntDir = NONE;
   seedlingID=0;
}


// -------------------------------------------------------------------------- //
// Assignment operator  - JM 06Dec2007                                        //
// To assign one object to another                                            //
// -------------------------------------------------------------------------- //

const Organism & Organism::operator=(const Organism & organismToCopy)
{
   if (this != &organismToCopy)
   {
      stock_ = organismToCopy.stock_ ;
      age = organismToCopy.age;
      ageRecruit = organismToCopy.ageRecruit;
      ageGameto = organismToCopy.ageGameto;
      ageMaturity = organismToCopy.ageMaturity;
      sporeStock = organismToCopy.sporeStock;
      species = organismToCopy.species;
      gametophyte = organismToCopy.gametophyte;
      gamMature = organismToCopy.gamMature;
      tempType = organismToCopy.tempType;
      shuntDir = organismToCopy.shuntDir;
   
      setNext(organismToCopy.getNext());
      setCell(organismToCopy.getCell());

   }
   return *this;
}


// -------------------------------------------------------------------------- //
// Initialiaze the static parameters pointer                                  //
// -------------------------------------------------------------------------- //

void Organism::initStaticPtr( const BaitParameters *inPtrParams, const EngineParameters * inPtrEngParams,
Fabric<Organism> *inPtrFabric )
{
   baitParamsPtr_ = inPtrParams ;
   engParamsPtr_ = inPtrEngParams;
   fabricPtr_     = inPtrFabric ;
}


//
// JM - 20Feb2007 Method to sample from Gaussian (Normal) distribution using GSL
//        (GNU Scientific Library).
//
double Organism::sampleGaussian(double stdDev)
{
   double sample=0.0;
   stdDev = fabs(stdDev);
   sample = gsl_ran_gaussian(Bait::rng, stdDev);

   // Correct for samples outside biologically realistic range.
   while (fabs(sample) > stdDev*4)
      sample = gsl_ran_gaussian(Bait::rng, stdDev);

   return sample;
}


//
// 30Oct2014 - calc mean gametophyte growth rate based on temp, day length, solar irradiance
//
void Organism::updGametGrowth(float currTemp, float dl, float solar)
{
   unsigned numLoopPerDay = 24;
            
   // 1. Base growth rate is function of water temperature (Morita, 2003):
   //    - @DL=12h, irrad=50 micromol m^-2 s^-1
   // Note, if stock starts at 0.05 & temp = 16C:
   //   - results in maturity after 14 days.
   //   - corresponds with 14 days from Choi et al., 2005 (DL=12,temp=16C,irradiance=60)
   // Thermal performance curve - Stevenson (1976) Eq (see Bulté & Blouin-Demers, 2006):
   // R^2=0.99983
   float K1 = 35.67618, 
         K2 = 0.15823, 
         K3 = 0.0147965,
         CTmin = 4.44542,
         CTmax = 28.23676,
         S = 10.62699;
   float A = 1 / (1 + K1 * exp(-K2 * (currTemp - CTmin)) );
   float B = 1 - exp( K3*(currTemp-CTmax) );   
   gametGrowth = (S * A) * B;

   // Correct for number of loops per day   
   gametGrowth = log(1.0+gametGrowth)/numLoopPerDay;   
   gametGrowth += sampleGaussian(gametGrowth*0.05);
 
}

//
// 30Oct2014 - calc mean sporophyte growth rate based on temp, day length, relative solar radiation
//           - see "Sporophyte Growth_v0.4.6_15Jan15.xlsx"
//
void Organism::updSporoGrowthMod(float currTemp, float dl, float solar)
{
   // reset counters for outputing relative growth rate
   sumSporoGrowthMod = 0.0;
   countSporoGrowthMod = 0;     
   sporoGrowthMod = 1.0;
   
   // 1. Day Length: 
   //    - Growth rate in response to day length (Pang & Luning, 2004, Fig. 7)
   //    - Meistematic zone is regulated by photoperiod
   // Hyperbolic equation of Jasby & Platt (1976):
   float P_maxDL = 1.5585;
   float alphaDL = 0.134154;
   float I_cDL = 0.0;      
   float actualDL = dl;
   
   sporoGrowthMod *= P_maxDL * ( 1.0 - exp(-alphaDL * (actualDL-I_cDL) / P_maxDL) );
      
   // 2. Water Temperature (Morita, 2003b (Fig.1)):
   //    - see also Saito, 1958(Fig.2); Gao, 2013(Fig.6)
   //    - Let temp modifier = 1.0 at 15 degrees celcius   
   // Thermal performance curve - Stevenson (1976) Eq (see Bulté & Blouin-Demers, 2006):
   // R^2=0.990421
   float K1 = 21.08749, 
         K2 = 0.212871, 
         K3 = 5.93443E-05,
         CTmin = 1.61626,
         CTmax = 28.27667,
         S = 3030.867;
   float A = 1 / (1 + K1 * exp(-K2 * (currTemp - CTmin)) );
   float B = 1 - exp( K3*(currTemp-CTmax) );   
   float thermal_fitness = (S * A) * B;
   
   sporoGrowthMod *= thermal_fitness;   
   sporoGrowthMod += sampleGaussian(sporoGrowthMod*0.05);
   
}

//
// 30Oct2014 - calc probability of gametophyte maturity based on temp, day length.
//           - see "Gametophyte growth_v0.5_16Jan15.xlsx"
//
void Organism::updProbGamMature(float currTemp, float dl)
{
   //unsigned numLoopPerDay = 24;
   float u=0.0, v=0.0, dlTempMod=0.0;
   float k=0.822, x0=17.61;

   // Probability of forming a new sporophyte
   probGamMature = 1.0;//baitParamsPtr_->probFertilise;     // TO DO: make proportional to num gametophytes in patch?

   // 1. Day length probability function    
   // Weibull curve fitted to data from Choi et al. (2005)
   //    - Table 1, average of all settlement densities   
   u = Organism::heightOfWeibull(dl, 10.96, 4.544) * 7.562;      
   
   // 2. Temperature effect
   // Sources: Choi (2005) Tab 1 and Morita (2003) Fig. 7      
   // Logistic curve
   v = 1.0 - ( 1.0 / (1.0 + exp(-k*(currTemp-x0))) );
   
   dlTempMod = u * v;        
   probGamMature *= dlTempMod;
   
   probGamMature += sampleGaussian(probGamMature*0.05);
}


// -------------------------------------------------------------------------- //
// Initialiaze the World (i.e. plate) static parameters pointer               //
// -------------------------------------------------------------------------- //

void Organism::initWorldPtr( World<Patch> * inWorld )
{
   ptrWorld = inWorld;
}


//
//   JM 13Jun2006 - Access Weibull distribution with specified parameters.
//      - see www.weibull.nl/weibullstatistics.htm for original equations
//
double Organism::heightOfWeibull(double x, double l_scale, double k_shape)
{
   double firstNumerator,
          firstDenominator,
          firstTerm,
          secondTerm;
   double yVal=0;

   firstNumerator = pow(x, k_shape-1);
   firstNumerator *= k_shape;
   firstDenominator = pow(l_scale, k_shape);
   firstTerm = firstNumerator / firstDenominator;
   
   secondTerm = pow((x/l_scale), k_shape);
   secondTerm *= -1;
   secondTerm = exp(secondTerm);
   yVal = firstTerm * secondTerm;
   return yVal;
}



// -------------------------------------------------------------------------- //
// Destructor                                                                 //
// -------------------------------------------------------------------------- //

Organism::~Organism()
{

}


// -------------------------------------------------------------------------- //
// Initialize all parameters                                                  //
// -------------------------------------------------------------------------- //

void Organism::init(int index, int rank)
{
   initAttributes();

   // initialize stock
   //    - when DL=12H, temp=16C, Irradiance = 60 micromol m^-2 s^-1
   stock_ = 0.05;
   stock_ += sampleGaussian(stock_*0.2);
   
   // total number of spores produced in season
   sporeStock = baitParamsPtr_->totalSpore; 
   sporeStock += sampleGaussian(sporeStock*0.2);  
   
   age = 0;
   ageRecruit = 0;
   ageGameto = 0;
   ageMaturity = 0;   
   gametophyte = true;     // starts off in gametophyte stage
   gamMature = false;
   releaseSpore = false;
   seedlingID = (rank * 1000) + index;
   numSeedlings++;
}

// -------------------------------------------------------------------------- //
// Initialize all parameters - supply stock as argument                       //
// -------------------------------------------------------------------------- //
void Organism::init(int index, int rank, double inStock)
{
   initAttributes();

   // initialize stock from input parameter
   stock_ = inStock + sampleGaussian(inStock*0.1);
   
   // total number of spores produced in season
   sporeStock = baitParamsPtr_->totalSpore; 
   sporeStock += sampleGaussian(sporeStock*0.2); 
    
   age = 0;
   ageRecruit = 0;
   ageGameto = 0;
   ageMaturity = 0;
   gametophyte = false;    // if stock is specified then assumed to be sporophyte
   gamMature = false;
   releaseSpore = false;
   seedlingID = (rank * 1000) + index;
   numSeedlings++;
}


// -------------------------------------------------------------------------- //
// The organism moves                                                         //
// -------------------------------------------------------------------------- //
void Organism::move()
{
   ( (Patch*) getCell() )->moveOrganism( this, NONE ) ;
}


//
// JM 23Jan14 - Growth algorithm - grow according to growthRate input parameter
//
void Organism::growth(unsigned nbLoops)
{
   unsigned numLoopPerDay = 24;
   
   // probability of early death proportional to size  
   if (getStock() >= 0.1)
   {
      if ( prematureDeath() )
         return;
   }

   if (gametophyte == true)
   {
      gametophyteGrowth();     
   }
   // Sporophyte growth (time to maturity, sporophyll production)
   else
   {
      float irradSurface = 575.0 * ((Plate*) ptrWorld)->getSolarRad();
      
      // irradiance decreases exponentially with depth
      // E(z) = E(0)e^-kz, where z=depth, k=attenuation coefficient
      float k = baitParamsPtr_->kdPAR;
      float irradiance = irradSurface * exp(-k*baitParamsPtr_->meanDepth);
      
      // plant length (microns) equals (stock * 3.0E+04)
      float plantL = getStock()*3.0E+05;
      
      // Sporophyte Base growth Rate:
      //    see "Sporophyte Growth_v0.4.6_15Jan15.xlsx"
      //    - fitted to data from:
      //    - Pang & Wu 1996, Fig. 6 (microscopic growth rate, 15deg,irradiance40)
      //    - Choi 2007, Table 2 (macroscopic growth rate, 15.6deg, irradiance unknown)
      //    Note: Stock 1.0 = 3.0E+05 microns in length
      float baseGrowth = 1.648 * pow(plantL,-0.242);
     
      
      // Photosynthesis rate is function of plant size:
      float P_max = 1.0;//0.514 * pow(plantL, 0.177);//0.453 * log(plantL) - 1.082;
      float alpha = 0.074 * pow(plantL, -0.183); //0.505 * pow(plantL, -0.342);
      float I_c = 1.4833 * pow(plantL, 0.1681);//1.94 * log(plantL) - 11.98;
      
      // Modify alpha because photosynthetic efficiency reduced in community
      // see (Binzer & Middelboe, 2005)
      alpha *= baitParamsPtr_->photo_alpha;     
      
      // Avoid extreme conditions when stock is very small
      // Use Choi (2005) data (on gametophytes) for growth of microscopic sporophytes
      if (P_max < 1.0)
         P_max = 1.0;
      if (alpha > 0.14224)
         alpha = 0.14224;
      if (I_c < 0.0)
         I_c = 0.0;
      
      float solarMod = P_max * ( 1.0 - exp(-alpha*(irradiance-I_c)/P_max) );
      
      //cout << plantL << ", " << solarMod << ", " << irradiance << endl; 
      float actualGrowth = baseGrowth * solarMod;  
      
      // modify according to other environmental conditions (temp, day length)
      actualGrowth *= sporoGrowthMod;

      // Track sporophyte growth rate in response to light, temp, day length
      sumSporoGrowthMod += (solarMod * sporoGrowthMod);
      countSporoGrowthMod++;
      if (countSporoGrowthMod > 0)
         avgSporoGrowthMod = sumSporoGrowthMod/(float)countSporoGrowthMod;
     
      actualGrowth += sampleGaussian(actualGrowth*0.05);
      
      if (actualGrowth > 0.7)
         actualGrowth = 0.7;        // max growth rate (from Shao-jun, 1996)
      if (actualGrowth < 0.0)
      {
         if (gsl_rng_uniform(Bait::rng) < 0.001)
         {
            (*this).deleteOrganism();    // Dies if growth rate is zero
            return;
         }
      }
      
      actualGrowth += 1.0;
      // convert from day^-1 to hour^-1
      actualGrowth = pow(actualGrowth, 1.0/numLoopPerDay);
           
      if (getSporeStock() > 1.0E-06)
         stock_ *= actualGrowth;      // exponential growth
      else
      {
         //resetAge();
         (*this).deleteOrganism();    // Dies if all spores released
         return;
      }
      
      // keep track of age (in loops) since first visible (used to id new recruits)
      age++;
      if (stock_ > baitParamsPtr_->macroSize)  
         ageRecruit++;
      
      // Record the age at maturity   
      if (stock_ > 1.0 && ageMaturity == 0)
      {
         ageMaturity = ageRecruit;
         sumAgeMaturity += ageMaturity;
         countMature++;
      }
   }
}

//
//
//
void Organism::gametophyteGrowth()
{
   ageGameto++;
   float irradSurface = 575.0 * ((Plate*) ptrWorld)->getSolarRad();
   float dl = ((Plate*) ptrWorld)->getCurrDayLength();

   // Modify growth rate according to irradiance
   // Hyperbolic function fitted to data from Choi et al., 2005
   // P_max (max level) and alpha (slope) function of day length
   // see "Photosynthesis vs irradiance curves_26Nov14.xlsx"
   float P_max = 0.29178 * exp(0.11224 * dl);
   float alpha = 2.849E-02 * dl - 0.1980;
   float I_c = 0.0;
   
   // irradiance decreases exponentially with depth
   // E(z) = E(0)e^-kz, where z=depth, k=attenuation coefficient
   float k = baitParamsPtr_->kdPAR;
   float irradiance = irradSurface * exp(-k*baitParamsPtr_->meanDepth);
   
   // Hyperbolic equation for photosynthesis vs. irradiance
   // see Campbell et al., 1999 and Jasby & Platt, 1976
   float newGrowth = P_max * ( 1.0 - exp(-alpha*(irradiance-I_c)/P_max) );
   newGrowth *= gametGrowth;   // peak rate, temp=20-25 (Morita et al., 2003)

   // Track gametophyte growth rate in response to light, temp, day length
   avgGametGrowthMod = newGrowth;
      
   // Add some random noise
   newGrowth += sampleGaussian(newGrowth*0.1);

   // 3. Overcrowding (sporophytes) inhibits growth
   // TO DO: competition for light
            
   // Apply new growth 
   if (newGrowth > 1.0E-06)
      stock_ += (stock_ * newGrowth);
   else
   {
      if (gsl_rng_uniform(Bait::rng) < 0.001)
      {
         (*this).deleteOrganism();    // Dies if growth rate is zero
         return;
      }
   }
             
   // Sporophyte is "born" under appropriate temperature conditions
   if (stock_ > 1.0)
   {
      // probMature: lower value means more gradual increase in rate of recruitment
      double probMature = 0.05;
      probMature *= probGamMature;
      
      if (gamMature == false)
      {
         //cout << "HERE: " << probGamMature << endl;
         if (gsl_rng_uniform(Bait::rng) < probMature)
         
            gamMature = true;
      }
      // chance that mature gametophyte is fertilised
      else
      {
         double probFert = 0.0;
         
         probFert = baitParamsPtr_->probFertilise * probGamMature;
         
         if (probFert > 0.0 && gsl_rng_uniform(Bait::rng) < probFert)
         {
            // Form sporophyte
            tryToReproduce();
            // Not able to produce new egg immediately (0.5 = arbitrary value)
            gamMature = false;
            stock_ = 0.5;
         }
      }
   }
}


//
// - Produce new sporophyte agent from gametophyte reproduction
//
void Organism::tryToReproduce()
{
   Organism   *ptrOrganism ;
   Direction_t dir;
   int rank = seedlingID / 1000;
   int i = seedlingID % 1000;

   // Create and initialise new sporophyte agent
   ptrOrganism = fabricPtr_->newAgent() ;
   ptrOrganism->init(i, rank, baitParamsPtr_->startStock);
   ptrOrganism->resetAge();
   ptrOrganism->setGametophyte(false);
 
   // Add to current patch
   dir = NONE;
   ( (Patch*) getCell()->getNeighbour(dir, true))->addOrganism( ptrOrganism );    
}



bool Organism::prematureDeath()
{
   double probEarlyDeathPow=0.0;
   float currDay = ((Plate*) ptrWorld)->getCurrDay();
   
   // 1. Probability of early death 
   // Type III Survivorship curve:
   //    - Weibull distribution for prob of early death:
   //       - deaths more concentrated among young agents
   //    - Data from Voisin thesis (2007):
   //       - Brest harbour 2004
   //       - see "Weibull_Voisin data_Nov14.xlsx"
   float ageInMonths=0;
   if (gametophyte == true)
      ageInMonths = (float)ageGameto/720.0;
   else
      ageInMonths = (float)ageRecruit/720.0;
    
   // If >6 months, switch to Type II (constant) survivorship curve
   if (ageInMonths > 6.0)
      ageInMonths = 6.0;           
   if (ageInMonths < 1.0)
      ageInMonths = 1.0;
   
   // Weibull function parameters   
   float k_shape = 0.16348;      // alpha
   float lambda_scale = 8.2E-06; // beta
   
   // Cosine function parameters (based on Brest harbour data):
   // see "Brest Survivorship Data_13Jul15.xlsx"
   float units = 0.0172142063;
   float k_amp = 0.2051320229;
   float k_trans_x = 248.5453109553;
   float k_trans_y = 2.6067436427;
   float l_amp = 0.2;
   float l_trans_x = 363.2100258048;
   float l_trans_y = 1.1654844907;
   
   // Derive k_shape cos function by transforming light cos function
   // Note: lambda remains the same for all locations (based on Brest data)
   // see "UGEN0.6.1_liget vs weibull_24Jul15.ods"
   k_amp = baitParamsPtr_->solarAmp * 0.02266447;
   k_trans_x = baitParamsPtr_->solar_c + 78.9546645;
   k_trans_y = (baitParamsPtr_->solar_d * 0.02266447) + 2.4016116;
 
   k_shape = k_amp * cos((currDay-k_trans_x) * units) + k_trans_y;
   lambda_scale = l_amp * cos((currDay - l_trans_x) * units) + l_trans_y;
   
   if (k_shape > exp(1))
      k_shape = exp(1);     
   
   if (gametophyte == true)
      probEarlyDeathPow = gsl_ran_weibull_pdf(ageInMonths, 8.2E-06, 0.164348) * 544.25;
   else
      probEarlyDeathPow = gsl_ran_weibull_pdf(ageInMonths, lambda_scale, k_shape);
   
   double propSurvive = 1.0 - probEarlyDeathPow;
   probEarlyDeathPow = (1.0 - pow(propSurvive,(1.0/720.0)));          
   
   if ( gsl_rng_uniform(Bait::rng) < probEarlyDeathPow)
   {
      (*this).deleteOrganism();
      return true;
   }
   else
      return false;
}


//
// JM 14Nov2013 - Produce spores and release into environment
//
void Organism::produceSpore()
{
   float currTemp = ((Plate*) ptrWorld)->getWaterTemp();
   unsigned numLoopPerDay = 24;
   
   // Only release spores if underwater and temperature above ~12-14 deg (Adance of Phycology in Japan, pg.307)
   if ( 
         gametophyte == false && releaseSpore == true 
         && getSporeStock() > 0.01 
         && ((Patch*)getCell())->getWaterDepth() > 0.0 )
   {
      // See Lewis et al, 1999 (PMID=10499277)
      float sporesReleased=0;

      sporesReleased = baitParamsPtr_->sporeProd;
      sporesReleased += sampleGaussian(sporesReleased*0.1);
      
      // To Do: Spore released is proportional to sporophyll length
      // Pang & Luning, 2004 
      
      // Saito pg 307, Advance of Phycology in Japan:
      // - Peak spore release at 17 - 22 degrees celcius
      // - Discharge complete after 20-40 days for medium sized sporophyll
      // Input parameter sporeRel is peak rate (i.e. at 17-22 degrees, i.e. 20 days)
      // Suto 1952: counts data for temperature versus spore shedding
      //    - Fitted logistic curve to data
      //    - see "Sporophyte Growth_v0.4.6_15Jan15.xlsx"
      float k=2.115, x0=13.519;
      float tempMod =  1.0 / ( 1.0 + exp(-k*(currTemp-x0)) );
     
      // Rate of release proportional to water temperature
      sporesReleased *= pow(tempMod,2);
      
      if(sporesReleased > 0.0)
      {
         if (sporesReleased < sporeStock)
         {
            ((Patch*)getCell())->addSpore(sporesReleased);        
            sporeStock -= sporesReleased;
         }
         else
         {
            ((Patch*)getCell())->addSpore(sporeStock);
            setSporeStock(0.0);
         }    
      } 
   }
   // Maturity when size is 1.09 (i.e. length>32.66cm, see Voisin (2007) Table 4.3)
   else if ( gametophyte == false && releaseSpore == false && stock_ >= 1.09)
   {
      // If conditions are correct will release spores (Suto, 1952)
      // Suto 1952: counts data for temperature versus spore shedding
      //    - Fitted logistic curve to data
      //    - see "Sporophyte Growth_v0.4.6_15Jan15.xlsx"
      float k=2.115, x0=13.519;
      float probRelease = 1.0 / ( 1.0 + exp(-k*(currTemp-x0)) );
      probRelease = 1.0 - pow((1.0-probRelease), (1.0/(10.0*numLoopPerDay)));
      
      if (gsl_rng_uniform(Bait::rng) < probRelease)
         releaseSpore = true;
   }
}




int Organism::round(double x)
{
    return static_cast<int>(floor(x + 0.5));
}

int Organism::round(float x)
{
    return static_cast<int>(floor(x + 0.5));
}


//
//   JM 13June - set general attributes of organism (species, temperature type) - not used
//
void Organism::initAttributes()
{

   species = U_PINNATIFIDA;         // Species
   tempType = MESO;            // Psychro-, Meso- or Thermo-phile. 
}


// -------------------------------------------------------------------------- //
// Determine the next patch to move in                                        //
// JM 13Jan2005 - implemented "run & tumble" model of bacterial movement   //
// -------------------------------------------------------------------------- //

Direction_t Organism::nextDirection()
{
   Direction_t   nextDir           ;

   // Retrieve next toward direction
   nextDir = getNextTowardDir() ;

   return nextDir ;
}



// -------------------------------------------------------------------------- //
// Apply shifts to organism position and retrieve new occupied patch          //
// -------------------------------------------------------------------------- //

Direction_t Organism::getNextTowardDir()
{
   Direction_t   nextDir;

   nextDir = (Direction_t) gsl_rng_uniform_int(Bait::rng,9);

   return nextDir ;
}


//
// JM 13Jun2006 - modifies intake rate according to pH and temperature of patch.
//
double Organism::phTempEffect(double intake)
{
   // [NOT USED]
   return intake;
}



//
// JM 13Jun2006 - Weibull distribution used to represent fitness effect of temperature.
//
double Organism::weibullDistribution(float xVal, TempType tempType)
{
   double peak, 
          percentage,
          returnVal=0.0, 
          tempModifier=0.0;   

   switch (tempType)
   {
      case PSYCHRO:
         peak = 0.17728;
         if(xVal > 0.0 && xVal < 20.0)
         {
            percentage = xVal/20;
            returnVal = (10.11 * (1-percentage))+.17;
         }
         tempModifier = Organism::heightOfWeibull(returnVal, 5, 2.1)/peak;
         break;
      case MESO:
         peak = 0.11498;
         if(xVal > 13.0 && xVal < 46.0)
         {
            percentage = (xVal-13)/33;
            returnVal = (15.0 * (1-percentage)) + .12;
         }
         tempModifier = Organism::heightOfWeibull(returnVal, 7, 1.8)/peak;
         break;
      case THERMO:
         peak = 0.1225;
         if(xVal > 42.0 && xVal < 76.0)
         {
            percentage = (xVal-42)/34;
            returnVal = (13.86 * (1-percentage))+.25;
         }
         tempModifier = Organism::heightOfWeibull(returnVal, 7, 2.0)/peak;
         break;

      case NONE_T:
         tempModifier = 1.0;
         break; 
   }
   return tempModifier;
}


//
// JM 12Mar2007 - Remove agent from patch
//
void Organism::deleteOrganism()
{
   // Only count mature sporophytes
   if (ageRecruit > 0 && ageMaturity >0)
   {
      sumAgeDeath += ageRecruit;
      countDead++;
   }

   // 1. Remove agent from patch
   ((Patch*) getCell())->removeOrgFromPatch(stock_);
  
   // 2. Set stock to zero (agent is dead)
   setStock(0.0);
   
   // 3. Reset spore stock to default value
   sporeStock = baitParamsPtr_->totalSpore; 
   sporeStock += sampleGaussian(sporeStock*0.2);
   
   resetAge();
   
   // Reset to default values:
   age = 0;
   ageRecruit = 0;
   ageMaturity = 0;   
   gametophyte = true;     // starts off in gametophyte stage
   gamMature = false;
   releaseSpore = false;

   // 4. Re-initialise agent/cell pointers to null
   Agent::setNext(NULL);
   Agent::setCell(NULL); 
}

//
// Set up a new Organism agent with stock, expSeed and seedlingID
//
void Organism::setNewTraits(double inStock, bool inRelSpore, int inSeedlingID)
{
   setStock(inStock);
   resetAge();
   setReleaseSpore((bool) inRelSpore);
   setSeedlingID(inSeedlingID); 
   
   // Default initial values:
   age = 0;
   ageRecruit = 0;
   ageMaturity = 0;   
   gametophyte = true;     // starts off in gametophyte stage
   gamMature = false;
   releaseSpore = false; 
}



// ------------------------------ End Of File ------------------------------- //
