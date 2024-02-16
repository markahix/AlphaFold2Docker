from flask import Flask, flash, request, redirect, render_template, session, Response

import uuid
CURRENT_JOBS={}

def random_job_identifier():
    id = uuid.uuid4()
    return id.hex