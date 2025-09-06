# Testing Guide for OBS-for-CRLB-ISAC Python Port

## Overview

This project uses **pytest** as the testing framework. This document explains the differences between unittest and pytest, and provides guidance on writing and running tests.

## unittest vs pytest

### unittest (Python's Built-in Framework)

**unittest** is Python's standard library testing framework:

```python
import unittest

class TestExample(unittest.TestCase):
    def setUp(self):
        self.config = SimulationConfig()
    
    def test_default_values(self):
        self.assertEqual(self.config.Nth, 4)
        self.assertTrue(self.config.Nt > 0)
        self.assertAlmostEqual(self.config.Pt, expected_value, places=5)

if __name__ == "__main__":
    unittest.main()
```

**Characteristics:**
- ✅ Built into Python (no extra installation)
- ❌ More verbose syntax
- ❌ Class-based approach required
- ❌ Less flexible fixtures
- ❌ Requires `self.assertEqual()`, `self.assertTrue()`, etc.

### pytest (Modern Testing Framework)

**pytest** is a third-party framework that's more Pythonic:

```python
import pytest

@pytest.fixture
def config():
    return SimulationConfig()

def test_default_values(config):
    assert config.Nth == 4
    assert config.Nt > 0
    assert abs(config.Pt - expected_value) < 1e-5

# Run with: pytest test_file.py
```

**Characteristics:**
- ✅ Simple, clean syntax
- ✅ Function-based approach (classes optional)
- ✅ Plain `assert` statements
- ✅ Powerful fixtures system
- ✅ Better parameterization
- ✅ Superior error reporting
- ✅ Rich plugin ecosystem

## Why We Use pytest

Our project uses pytest because:

1. **Cleaner Syntax**: No need for `self.assertEqual()` - just use `assert`
2. **Better Fixtures**: More flexible setup/teardown with `@pytest.fixture`
3. **Parameterization**: Easy to test multiple scenarios with `@pytest.mark.parametrize`
4. **Better Output**: More informative error messages and test reports
5. **Integration**: Works well with existing Python tools and CI/CD

## Test Structure

Our tests are organized as follows:

```
python-port/tests/
├── calculate_fim/
│   └── test_calculate_fim.py
├── construct_matrixQ/
│   └── test_construct_matrixQ.py
├── steering_matrix/
│   └── test_validation.py
├── test_compare_mat_files.py
├── test_simulation_config.py
└── TESTING_GUIDE.md
```

## Running Tests

### Prerequisites

Follow [python-port/README.md](../README.md) for virtual environment, dependencies setup.

### Running Tests

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest test_simulation_config.py

# Run specific test class
pytest test_simulation_config.py::TestSimulationConfig

# Run specific test method
pytest test_simulation_config.py::TestSimulationConfig::test_default_configuration

# Run tests with coverage
pytest --cov=utils

# Run tests matching a pattern
pytest -k "config"
```

## Writing New Tests

### Basic Test Function

```python
def test_my_feature():
    """Test description"""
    # Arrange
    config = SimulationConfig(Nth=8)
    
    # Act
    result = config.some_method()
    
    # Assert
    assert result == expected_value
```

### Using Fixtures

```python
@pytest.fixture
def sample_config():
    """Fixture providing a test configuration"""
    return SimulationConfig(
        Nth=4, Ntv=4, K=2, M=1,
        random_seed=42
    )

def test_with_fixture(sample_config):
    assert sample_config.Nt == 16
```

### Parameterized Tests

```python
@pytest.mark.parametrize("input_val,expected", [
    (0, 1.0),
    (10, 10.0),
    (-10, 0.1)
])
def test_db2pow(input_val, expected):
    assert abs(db2pow(input_val) - expected) < 1e-10
```

### Testing Exceptions

```python
def test_invalid_input():
    with pytest.raises(ValueError, match="Invalid parameter"):
        SimulationConfig(Nth=-1)
```

### Testing Output

```python
def test_print_output(capsys):
    demo_function()
    captured = capsys.readouterr()
    assert "Expected output" in captured.out
```

## Best Practices

1. **Test Names**: Use descriptive names starting with `test_`
2. **Docstrings**: Add brief descriptions for complex tests
3. **Arrange-Act-Assert**: Structure tests clearly
4. **Isolation**: Each test should be independent
5. **Edge Cases**: Test boundary conditions and error cases
6. **Fixtures**: Use fixtures for complex setup
7. **Mocking**: Mock external dependencies

## Example: Converting unittest to pytest

### Before (unittest):
```python
import unittest

class TestConfig(unittest.TestCase):
    def setUp(self):
        self.config = SimulationConfig()
    
    def test_antennas(self):
        self.assertEqual(self.config.Nth, 4)
        self.assertEqual(self.config.Nt, 16)
```

### After (pytest):
```python
import pytest

@pytest.fixture
def config():
    return SimulationConfig()

def test_antennas(config):
    assert config.Nth == 4
    assert config.Nt == 16
```

## Integration with Project

The tests integrate seamlessly with the existing project structure:

- **42 total tests** (30 new + 12 existing)
- **All tests pass** ✅
- **Consistent with existing pytest conventions**
- **Comprehensive coverage** of the SimulationConfig system
- **Easy to extend** for future modules

## Continuous Integration

For CI/CD pipelines, use:

```yaml
# .github/workflows/tests.yml
- name: Run tests
  run: |
    python -m virtual-environment venv
    source venv/bin/activate
    pytest -v --cov=utils --cov-report=xml
```

This ensures all tests pass before merging code changes.
