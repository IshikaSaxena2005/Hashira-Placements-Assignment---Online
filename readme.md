# Hashira Placements Assignment - Polynomial Solver

This project solves polynomial coefficient problems using Lagrange interpolation, specifically designed for the Hashira Placements Assignment.

## Problem Description

Given a set of points with coordinates in different number bases, the task is to:
1. Convert the values from their respective bases to decimal
2. Use Lagrange interpolation to reconstruct the polynomial
3. Find the constant term (coefficient of x^0) of the polynomial

## Project Structure

```
hashira-placements-solver/
├── polynomial_solver.js      # Main solver implementation
├── package.json             # Node.js package configuration
├── create_test_files.js     # Script to create test case JSON files
├── testcase1.json          # First test case (created by script)
├── testcase2.json          # Second test case (created by script)
└── README.md               # This file
```

## Setup and Installation

1. Make sure you have Node.js installed on your system
2. Create a new directory for the project:
   ```bash
   mkdir hashira-placements-solver
   cd hashira-placements-solver
   ```
3. Copy all the provided files into this directory
4. Create the test case JSON files:
   ```bash
   node create_test_files.js
   ```

## Usage

### Running the solver with embedded test cases:
```bash
node polynomial_solver.js
```

### Running with external JSON files:
```javascript
const { solveFromFile } = require('./polynomial_solver');

// Solve from file
solveFromFile('testcase1.json');
solveFromFile('testcase2.json');
```

## Algorithm Explanation

The solution uses **Lagrange Interpolation** to reconstruct the polynomial:

1. **Base Conversion**: Convert all values from their respective bases to decimal
2. **Point Extraction**: Extract (x, y) coordinates where x is the key and y is the converted value
3. **Lagrange Interpolation**: Use the first k points to reconstruct the polynomial
4. **Constant Term**: Evaluate the polynomial at x=0 to get the constant term

### Lagrange Interpolation Formula
For a polynomial P(x) passing through points (x₀, y₀), (x₁, y₁), ..., (xₙ₋₁, yₙ₋₁):

```
P(x) = Σ(i=0 to n-1) yi * Li(x)

where Li(x) = Π(j=0 to n-1, j≠i) (x - xj) / (xi - xj)
```

To find the constant term, we evaluate P(0).

## Expected Output

For the given test cases:
- **Test Case 1**: The constant term should be calculated and displayed
- **Test Case 2**: The constant term should be calculated and displayed

## File Formats

### Input JSON Format:
```json
{
    "keys": {
        "n": 4,        // Number of roots provided
        "k": 3         // Minimum roots required (degree + 1)
    },
    "1": {
        "base": "10",
        "value": "4"
    },
    // ... more points
}
```

### Key Information:
- `n`: Total number of points available
- `k`: Minimum number of points needed (polynomial degree + 1)
- Each numbered key represents a point where the key is the x-coordinate
- `base`: The number base of the value
- `value`: The y-coordinate in the specified base

## Dependencies

- Node.js (built-in modules only: `fs`)
- No external npm packages required

## Notes

- The solution handles various number bases (binary, octal, decimal, hexadecimal, etc.)
- Uses JavaScript's built-in `parseInt()` function for base conversion
- Implements mathematical precision handling for large numbers
- Works with the Shamir's Secret Sharing cryptographic scheme principles

## Testing

The code includes both embedded test cases and the ability to read from external JSON files. Run the main script to see results for both test cases.