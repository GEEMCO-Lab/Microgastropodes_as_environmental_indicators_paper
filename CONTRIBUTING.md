# Contributing to Microgastropods as Environmental Indicators

Thank you for your interest in contributing to this research project! This document provides guidelines for contributing to the repository.

## How to Contribute

### Reporting Issues

If you find a bug, have a question, or want to suggest an improvement:

1. Check if the issue already exists in the [Issues](https://github.com/GEEMCO-Lab/Microgastropodes_as_environmental_indicators_paper/issues) section
2. If not, create a new issue with a clear title and description
3. Include relevant details such as:
   - What you expected to happen
   - What actually happened
   - Steps to reproduce the issue
   - Your R version and package versions

### Suggesting Enhancements

We welcome suggestions for improvements:

1. Open an issue describing the enhancement
2. Explain why this enhancement would be useful
3. Provide examples if possible

### Code Contributions

If you'd like to contribute code:

1. **Fork the repository**
   ```bash
   git clone https://github.com/GEEMCO-Lab/Microgastropodes_as_environmental_indicators_paper.git
   ```

2. **Create a new branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make your changes**
   - Follow the existing code style
   - Add comments where necessary
   - Update documentation if needed

4. **Test your changes**
   - Ensure all scripts run without errors
   - Verify that outputs are as expected

5. **Commit your changes**
   ```bash
   git commit -m "Brief description of your changes"
   ```

6. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

7. **Submit a pull request**
   - Provide a clear description of the changes
   - Reference any related issues

## Code Style Guidelines

### R Code Style

- Use `<-` for assignment, not `=`
- Use meaningful variable names
- Indent with 2 spaces
- Maximum line length of 80 characters when possible
- Add comments for complex operations
- Use `tidyverse` style guide: https://style.tidyverse.org/

Example:
```r
# Good
species_diversity <- calculate_diversity(community_matrix)

# Avoid
sd = calc_div(cm)
```

### File Organization

- Keep related functions together
- Use descriptive file names
- Add README files to explain directory contents
- Document any non-obvious design decisions

## Data Contributions

If you have data to contribute:

1. Ensure you have permission to share the data
2. Document the data collection methodology
3. Provide metadata (column descriptions, units, etc.)
4. Use standard formats (CSV for tabular data)
5. Include data provenance information

## Documentation

- Update README.md if you add new features
- Add comments to explain complex code
- Update the manuscript if results change
- Keep the requirements.txt and renv.lock files current

## Testing

Before submitting:

1. Run all scripts from a clean environment
2. Verify that all outputs are generated correctly
3. Check that the manuscript compiles without errors
4. Ensure no sensitive or private data is included

## Questions?

If you have questions about contributing, please:
- Open an issue with the "question" label
- Contact the project maintainers

## Recognition

All contributors will be acknowledged in the project documentation and, where appropriate, in the manuscript acknowledgments.

Thank you for helping make this research more reproducible and transparent!
