# Contributing to the Package development

## Contributing to development of the *ibis.iSDM* R-package

We welcome contributions to the *ibis.iSDM* R-package. These
contributions could be simple typo fixes, additions to the documentation
and **testthat** tests, enhancing the vignettes to provide a greater
understanding this package, or completely new functions.

For the latter, please get in touch with the package author or one of
the maintainers first. Pull requests to the *master* branch require a
confirmation and code review from package maintainers.

### Development guidelines

- The *ibis.iSDM* contains primarily functions for fitting models.
- Speed and flexibility are key
- Don’t repeat yourself. Create new functions and if necessary classes.
  Equally try to reuse common names from R, e.g. *plot*, *summary*
- Please run code *checks* and *tests* regularly.
- Avoid using additional package dependencies where possible.
- Comment your code!!
- Use assertions to verify inputs to functions.
- If bored, please write unit tests and ensure they all evaluate
  (CRTL+SHIFT+T)!

(also see [issues](https://github.com/iiasa/ibis.iSDM/issues) and
[projects](https://github.com/iiasa/ibis.iSDM/projects)) for open issues
