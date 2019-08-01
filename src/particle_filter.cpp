/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <limits>

#include "helper_functions.h"

#define EPS 0.00001 // Just a small number

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  //std::cout<<"[debug]here in init \n";

  // create random particles around init (x,y) with theta variance
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  num_particles = 100;  // TODO: Set the number of particles
  for(int i =0;i<num_particles;++i) {
      weights.push_back(1.0); // init all weights to 1
      
      // sample around the inital (x,y) to init particles
      Particle p;
      p.id = i;
      p.x = dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
      p.weight = 1.0;

      particles.push_back(p); // add to particles of all particles
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    //std::cout<<"[debug]here in the prediction \n";
   // loop through all particles and predict with noise from their
   // previous location to the new location
    std::default_random_engine gen;
    // add noise to the moving process with mean of the update position and angle
    // with standard deviation of the moving process provided
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

    for(auto& p: particles) {
        // compute the predited next position x , y and theta
        if(fabs(yaw_rate) < EPS) { // different precition formula depend if yaw_rate changes
            p.x += velocity*delta_t*cos(p.theta);
            p.y += velocity*delta_t*sin(p.theta);
        } else {
            p.x += velocity/yaw_rate*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
            p.y += velocity/yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate*delta_t));
            p.theta += yaw_rate*delta_t;
        }

        // add noise to the x, y and theta
        p.x += dist_x(gen);
        p.y += dist_y(gen);
        p.theta += dist_theta(gen);

    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  // predicted is the landmark that predicted from each particles
  // observations is the actual sernsor data of the landmark
  // associate the nearest  
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   // for each observeation in map coordinates find its nearest match
   // of landmark and store associated landmark index in the particle
   for(auto& ob: observations) {
       double small_dist = std::numeric_limits<double>::max(); // record smallest distantance
       int index = -1;
       for(auto& lm: predicted) {
           auto distance = dist(ob.x,ob.y,lm.x,lm.y);//sqrt(pow(ob.x-lm.x, 2) + pow(ob.y - lm.y, 2));
           if(distance < small_dist){
               small_dist = distance;
               index = lm.id;
           }
       }

       // record the nearest landmark index to the particle
       ob.id = index; // record the id of the landmark this observaton binds to
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

    // given oberservation in car's coordinates translate into the map coordinates
    // then do the associateion with the closet landmark for each transfromed observation
    // loop through each particle and transform obserations from its perspetive to
    // the map coordinates and associate each observation to the nearest map landmark

    std::default_random_engine gen;
    for(auto& p: particles) {

        // do the car coordinates to map coordinates transformation
        vector<LandmarkObs> observe;
        for(auto& ob: observations) {
           // transform observed data in car coordinates for each partical
            // into map coordinates
            auto xm = p.x + cos(p.theta)*ob.x - sin(p.theta)*ob.y;
            auto ym = p.y + sin(p.theta)*ob.x + cos(p.theta)*ob.y;

            // add to the observation vector
            LandmarkObs lob;
            lob.x = xm;
            lob.y = ym;
            observe.push_back(lob);
        }

        // find all landmarks that is within the landmark range for each particle
        vector<LandmarkObs> predicted;  
        for(auto& lm: map_landmarks.landmark_list){
            auto distance = dist(p.x, p.y, lm.x_f, lm.y_f); //sqrt(pow(p.x-lm.x_f, 2) + pow(p.y - lm.y_f, 2));
            if(distance <= sensor_range){ // only add map landmarks within sensor range of this particle
                LandmarkObs lob;
                lob.x = lm.x_f;
                lob.y = lm.y_f;
                lob.id = lm.id_i;

                predicted.push_back(lob);
            } 
        }
        
        // now we have predicted for each particle all landmarks within sensor range
        // observed is the transformed map cooridnates sensing data of landamarks
        // make each observation associate with the nearest landmark for this particle
        dataAssociation(predicted, observe);
        
        // now we have the nearest land mark index stored in the particle
        // same for the sensed data transformed in map coordinates
        // assign weight for this particle
        p.weight= 1;
        vector<double> sense_x, sense_y;
        vector<int> associations;
        for(auto& ob: observe){
            // find the associated id of land mark of this observation
            // it points to the landamarks that this observation assocaiated with
            auto it = std::find_if(predicted.begin(), predicted.end(),
                    [ob](LandmarkObs x){return x.id == ob.id;});
            
            auto one_weight= multiv_prob(std_landmark[0], std_landmark[1], 
                                         ob.x, ob.y,it->x,it->y);

            p.weight *= (one_weight < EPS) ? EPS : one_weight;

            // for debuging associate particles with observerations
            sense_x.push_back(ob.x);
            sense_y.push_back(ob.y);
            associations.push_back(ob.id);
        }
        SetAssociations(p, associations, sense_x, sense_y);
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    // get all weights 
    vector<double> new_weights;
    for(auto& p: particles) 
        new_weights.push_back(p.weight);

    std::vector<Particle> new_particles;
    // replace with prob according to weight
    std::default_random_engine gen;
    std::discrete_distribution<> dist(new_weights.begin(),new_weights.end());

    // replace with number of paricles times and update the particle
    for(int i = 0; i < num_particles;++i) 
        new_particles.push_back(particles[dist(gen)]);

    particles = new_particles; // update to the resampled particle
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {

  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
