Training a Topological Descriptor Model
========================================

In this document, we will walk through the steps to train a surrogate energy model for nanoparticles using Topological Descriptors.

## Tutorial Overview
This tutorial covers the following steps:

1. **Problem Introduction**
   - Define objectives and expected outcomes.

2. **Nanoparticle Creation**
   - Create nanoparticles using the `Nanoparticle` class and generate a training set.

3. **Feature Extraction**
   - Use the `TopologicalFeatureClassifier` class to extract features.

4. **Model Training**
   - Train a Bayesian Ridge Regression model with `BayesianRRCalculator`.

5. **Model Evaluation**
   - Evaluate model performance and visualize results.

6. **Model Saving**
   - Save the trained model for future use.

## Example Code

```python
# Import necessary modules from NPlib
from npl.core import Nanoparticle
from npl.calculators import EMTCalculator
from npl.descriptors.global_feature_classifier import testTopologicalFeatureClassifier
from npl.calculators import BayesianRRCalculator
from npl.utils.utils import plot_learning_curves
import numpy as np
import matplotlib.pyplot as plt
import pickle

# Function to create a training set of nanoparticles
def create_octahedron_training_set(n_particles, height, trunc, stoichiometry):
    emt_calculator = EMTCalculator(fmax=0.2, steps=1000)
    
    training_set = []
    for i in range(n_particles):
        p = Nanoparticle()
        p.truncated_octahedron(height, trunc, stoichiometry)
        emt_calculator.compute_energy(p)
        training_set.append(p)
        
    return training_set

# Create the training set with 40 particles
stoichiometry = {'Pt': 55, 'Au': 24}
training_set = create_octahedron_training_set(40, 5, 1, stoichiometry)

# Extract features
classifier = testTopologicalFeatureClassifier(list(stoichiometry.keys()))
for p in training_set:
    classifier.compute_feature_vector(p)

# Train the Bayesian Ridge Regression model
calculator = BayesianRRCalculator(classifier.get_feature_key())
calculator.fit(training_set, 'EMT', validation_set=0.1)

# Evaluate the model
X = [p.get_feature_vector('TFC') for p in training_set]
y = [p.get_energy('EMT') for p in training_set]
n_atoms = training_set[0].get_n_atoms()
plot_learning_curves(X, y, n_atoms, calculator.ridge, n_splits=10, train_sizes=range(4, 30, 2), y_lim=(0, 2))

# Visualize the coefficients
coefficients = calculator.get_coefficients()
feature_names = classifier.get_feature_labels()
plt.figure(figsize=(10, 6))
plt.bar(range(len(coefficients)), coefficients)
plt.hlines(0, 0, len(coefficients), linestyles='dashed')
plt.xticks(range(len(coefficients)), feature_names, rotation=90)
plt.xlabel('Coefficient Index')
plt.ylabel('Coefficient Value')
plt.title('Fitting Coefficients')
plt.show()

# Save the trained model
calculator.save('bayesian_rr_calculator.pkl')
