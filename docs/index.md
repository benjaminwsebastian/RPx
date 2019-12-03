---
title: "Getting started with RPx"
keywords: sample homepage
tags: [getting_started]
sidebar: mydoc_sidebar
permalink: index.html
summary: These brief instructions will help you get started with the RPx. The other topics in this help provide additional information and detail about working with other aspects of RPx.
---

{% include note.html content="This project is still in early development so it is not registered in the PyPi package index. It will be added upon the v1.0 release" %}

## Installation - dev

Follow these instructions to build the theme.

### 1. Install dependencies

To use RPx you need Python 3, Pip, Virtualenv, and Pytest

Virtualenv can be installed with ```pip install virtualenv```

Pytest can be installed with ```pip install -U pytest```

### 2. Clone repo

First, download or clone RPx from the [Github repo](https://github.com/benjaminwsebastian/RPx):

```git clone git@github.com:benjaminwsebastian/RPx.git && cd RPx``` (use git clone ```https://github.com/benjaminwsebastian/RPx && cd RPx```)

### 3. Create dev environment

Create a virtual environment:

```python3 -m venv venv```  

```source venv/bin/activate```  

```pip install -r requirements.txt```

### 4. Test installation

You can test the installation with ```pytest -v```

    Note that it should fail on one test

## Quick Start

{% include links.html %}
